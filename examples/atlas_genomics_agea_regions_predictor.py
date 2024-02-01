"""
In this example we train a neural network to predict the cosmos regions from the gene expression
This demonstrates how to align and sample brain regions relative to the gene expression volumes.

This examples requires `sklearn` and `seaborn` to be installed on top of the `iblatlas` requirements.
"""

import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import sklearn.metrics

from iblatlas.genomics import agea

# %% Load the agea atlas
df_genes, gene_expression_volumes, atlas_agea = agea.load()

# %% remap the the agea atlas at the cosmos level parcellation
ne = gene_expression_volumes.shape[0]
sel = atlas_agea.label.flatten() != 0  # remove void voxels
# reshape in a big array nexp x nvoxels this takes a little while
gexps = gene_expression_volumes.reshape((ne, -1))[:, sel].astype(np.float32).transpose()
aids = atlas_agea.regions.id[atlas_agea.label.flatten()[sel]]
aids_cosmos = atlas_agea.regions.remap(aids, 'Allen-lr', 'Cosmos')

# %% now we learn to predict the cosmos labels from the gene expression
X_train, X_test, y_train, y_test = train_test_split(gexps, aids)
scaler = StandardScaler()
scaler.fit(gexps)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
clf = MLPClassifier(random_state=1, max_iter=300, verbose=True).fit(X_train, y_train)
clf.predict_proba(X_test[:1])
clf.predict(X_test)
clf.score(X_test, y_test)
classes_labels = atlas_agea.regions.id2acronym(clf.classes_)

# %% Plot the confusion matrix
cm = sklearn.metrics.confusion_matrix(y_test, clf.predict(X_test), normalize='pred')
sklearn.metrics.ConfusionMatrixDisplay(cm, display_labels=classes_labels).plot(ax=plt.gca(), cmap='magma')

fig, ax = plt.subplots(1, 1, figsize=(7, 6))
sns.heatmap(cm.T * 100, vmin=0, vmax=10, cmap='Blues', annot=True, ax=ax, fmt='.1f')
ax.set(
    xticklabels=classes_labels,
    yticklabels=classes_labels,
    xlabel='True region',
    ylabel='Predicted region',
    title='Confusion Matrix (%)'
)
