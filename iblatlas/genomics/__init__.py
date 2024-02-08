"""A package for working with Allen genomics datasets: AGEA and MERFISH.

AGEA
----

This package provides a way to load the Allen Gene Expression volumes.
The 4345 volumes have been registered and formatted into a binary file.

     agea/
     ├── gene-expression.bin
     ├── gene-expression.pqt
     ├── image.npy
     └── label.npy

-   gene-expression.bin is a float-16 binary file containing the gene expression volumes.
In c-order, the dimensions are (4345, 58, 41, 67) that corresponds to (nexperiments, ml, dv, ap) at 200 um.
-   gene-expression.pqt is a parquet file describing the 4345 genes corresponding to the
gene expression volumes.
-   image.npy: the Allen atlas diffusion imaging volume downsampled at the gene expression resolution
-   label.npy: the Allen atlas region label volume downsampled at the gene expression resolution
See the building scripts in ./genomics/gene_expression_scrapping/05-generate-atlas-templates.py

[1] E. S. Lein et al., “Genome-wide atlas of gene expression in the adult mouse brain,”
 Nature, vol. 445, no. 7124, Art. no. 7124, Jan. 2007, doi: 10.1038/nature05453.
[2] L. Ng et al., “An anatomic gene expression atlas of the adult mouse brain,”
 Nat Neurosci, vol. 12, no. 3, Art. no. 3, Mar. 2009, doi: 10.1038/nn.2281.


MERFISH
-------

Spatially registered cell types from single cell transcriptomics data.

This package provides a way to load the MERFISH data from the Allen Brain Cell Atlas.
We formatted the original CSV files from the 2023-12-15 release into parquet files for faster loading and smaller hard
drive footprint.

    merfish/
    ├── genes.pqt
    ├── neurotransmitters.pqt
    ├── classes.pqt
    ├── subclasses.pqt
    ├── supertypes.pqt
    ├── clusters.pqt
    ├── C57BL6J-638850_cells.pqt
    ├── Zhuang-ABCA-1_cells.pqt
    ├── Zhuang-ABCA-2_cells.pqt
    ├── Zhuang-ABCA-3_cells.pqt
    └── Zhuang-ABCA-4_cells.pqt

-   *_cells.pqt: Each dataframe corresponds to a given subject. The concatenation of those 5 dataframes lead to
8_879_868, 11 cells with the following columns:
    -   'brain_section_label': the label of the brain section (subject and section): "Zhuang-ABCA-1.085"
    -   'donor_label': the label of the subject
     -  'neurotransmitter': neurotransmitter label {<NA>, 'Glut', 'Chol', 'GABA-Glyc', 'GABA','Dopa',
     'Glut-GABA', 'Hist', 'Sero', 'Nora'}
     -  'class': direct index of the class record in df_classes
     -  'subclass': direct index of the subclass record in df_subclasses
     -  'supertype': direct index of the supertype record in df_supertypes
     -  'cluster': direct index of the cluster record in df_clusters
     -  'x', 'y', 'z': coordinates of the cell in IBL space (see: iblatlas.atlas.AllenAtlas)
     -  'Allen_id': allen region unique identifier

The cells are classified hierarchically, from high level to low level: classes, subclasses, supertypes and clusters.
-   df_classes: a dataframe of classes (35, 3), where each record corresponds to a single class
-   df_subclasses: a dataframe of subclasses (339, 4), where each record corresponds to a single subclass
-   df_supertypes: a dataframe of supertypes (1202, 4), where each record corresponds to a single supertype
-   df_clusters: a dataframe of clusters (5323, 5), where each record corresponds to a single cluster

Additional metadata:
-   df_neurotransmitters: a dataframe of neurotransmitters (9, 2), index is the neurotransmitter label
-   df_genes: a dataframe of genes (1122), this could be used in conjunction with raw gene expressions data (not implemented)

[1] Z. Yao et al., “A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain,”
 Nature, vol. 624, no. 7991, Art. no. 7991, Dec. 2023, doi: 10.1038/s41586-023-06812-z.
[2] M. Zhang et al., “Molecularly defined and spatially resolved cell atlas of the whole mouse brain,”
 Nature, vol. 624, no. 7991, Art. no. 7991, Dec. 2023, doi: 10.1038/s41586-023-06808-9.

"""