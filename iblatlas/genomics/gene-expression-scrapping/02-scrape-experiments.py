from pathlib import Path
import subprocess
import shutil

import numpy as np
import pandas as pd
from zipfile import ZipFile
import skimage.io as io

folder_download = Path("/datadisk/gdrive/2022/08_gene_expression_atlas")
file_parquet = folder_download.joinpath('gene-expression.pqt')
folder_experiments = folder_download.joinpath('experiments')
experiments = pd.read_parquet(file_parquet)


# %% scrape all zip files
folder_experiments.mkdir(exist_ok=True)
for i, exp in experiments.iterrows():
    file_experiment = folder_experiments.joinpath(f"{exp.id}_{exp.gene}.zip")
    if not file_experiment.exists():
        print(i, file_experiment)
        cmd = f"wget http://api.brain-map.org/grid_data/download/{exp.id} -P {folder_experiments}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        info, error = process.communicate()
        shutil.move(folder_experiments.joinpath(exp.id), file_experiment)


# %% Interpret
file_out_npy = folder_download.joinpath("gene-expression.bin")
scratch_dir = folder_download.joinpath('scratch')

with open(file_out_npy, 'wb+') as fp:
    for i, exp in experiments.iterrows():
        print(i, file_experiment)
        file_experiment = folder_experiments.joinpath(f"{exp.id}_{exp.gene}.zip")
        scratch_dir.mkdir(exist_ok=True)
        with ZipFile(file_experiment, 'r') as zip:
            zip.extractall(scratch_dir)

        #  (58, 41, 67) corresponds to ML, DV, AP
        mhd_file = next(scratch_dir.glob('*.mhd'))
        img = io.imread(mhd_file, plugin='simpleitk')

        shutil.rmtree(scratch_dir)

        img.astype(np.float16).tofile(fp)
        # import matplotlib.pyplot as plt
        # plt.imshow(np.squeeze(img[:, :, 25].T))
