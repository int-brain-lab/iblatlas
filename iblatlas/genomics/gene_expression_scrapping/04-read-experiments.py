from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from ibllib.atlas import AllenAtlas, BrainCoordinates

DIM_EXP = (4345, 58, 41, 67)
ne, nml, ndv, nap = DIM_EXP
folder_download = Path("/datadisk/gdrive/2022/08_gene_expression_atlas")
file_parquet = folder_download.joinpath('gene-expression.pqt')
file_bin = file_parquet.with_suffix(".bin")

folder_experiments = folder_download.joinpath('experiments')
experiments = pd.read_parquet(file_parquet)

ba = AllenAtlas()

gexp_all = np.memmap(file_bin, dtype=np.float16, mode='r', offset=0, shape=DIM_EXP)

# %% (58, 41, 67) corresponds to ML, DV, AP
ie = np.where(experiments['gene'] == 'Pantr1')[0][0]

gexp = gexp_all[ie, :, :, :].astype(np.float32)
gexp[gexp <= 0] = np.nan
bc = BrainCoordinates(nxyz=[nml, nap, ndv], dxyz=[200, 200, 200])

xyz = np.array([0, 0, -0.002])
ba.plot_slices(xyz, volume="annotation")

iml, iap, idv = np.round(ba.xyz2ccf(xyz) / 200).astype(np.int32)

# BrainAtlas(image=gexp, label, dxyz, regions, iorigin=[0, 0, 0],
#              dims2xyz=[0, 1, 2], xyz2dims=[0, 1, 2]
fig, axs = plt.subplots(2, 2)
axs[0, 0].imshow(np.squeeze(gexp[:, :, iap]).T)
axs[0, 1].imshow(np.squeeze(gexp[iml, :, :]))
axs[1, 0].imshow(np.squeeze(gexp[:, idv, :]))


# %% This cell takes
ccf_coords = np.meshgrid(*[np.arange(DIM_EXP[i]) * 200 for i in [1, 2, 3]])  # ml, dv, ap
xyzs = ba.ccf2xyz(np.c_[ccf_coords[0].flatten(), ccf_coords[2].flatten(), ccf_coords[1].flatten()])
aids = ba.get_labels(xyzs, mode='clip', mapping='Cosmos')

sel = ~np.logical_or(aids == 0, aids == 997)  # remove void and root voxels
# reshape in a big array nexp x nvoxels
gexps = gexp_all.reshape((ne, nml * nap * ndv))[:, sel].astype(np.float32).transpose()
xyzs = xyzs[sel]
aids = aids[sel]

# %% ok ready for the learning part
X_train, X_test, y_train, y_test = train_test_split(gexps, aids, stratify=aids)
scaler = StandardScaler()
scaler.fit(gexps)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
clf = MLPClassifier(random_state=1, max_iter=300, verbose=True).fit(X_train, y_train)
clf.predict_proba(X_test[:1])
clf.predict(X_test)
clf.score(X_test, y_test)
