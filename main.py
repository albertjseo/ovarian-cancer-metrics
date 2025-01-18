import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt

# Set default matplotlib params
sc.settings.verbosity = 3  # providing additional information / feedback from program. 0 - errors, 1 - warning messages, 2 - informational message, 3 - hints, 4 - debugging
sc.logging.print_header() #Versions that might influence the numerical results. Matplotlib and Seaborn are excluded from this.
sc.settings.set_figure_params(dpi=100, facecolor="white") # sets resolutions, sizing, styling, and formatting of figure.

# adata is an AnnData object that can be sliced like a dataframe. AnnData stores a data matrix with observations, variables, and unstructures annotations.
# read_10x_mtx returns an AnnData object
adata = sc.read_10x_mtx(
    "/Users/albertseo/Development/ovarian-cancer-metrics/filtered_feature_bc_matrix",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names
    cache=True,  # write a cache file for faster subsequent reading
)

adata.var_names_make_unique() # apply common gene name

# Preprocessing
sc.pl.highest_expr_genes(adata, n_top=25) # plotting top 25 expressing genes from adata
sc.pp.filter_cells(adata, min_genes=200) # filter out cells with low gene counts
sc.pp.filter_genes(adata, min_cells=3) # basic quality control filter to remove genes expressed in very few cells

#filtered out 597 cells that have less than 200 genes expressed & 728 genes that are detected in less than 3 cells

# number of genes with positive counts in cell vs total number of cells
adata.var["mito"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], log1p=False, inplace=True)

sns.jointplot(
    data=adata.obs,
    y="total_counts",
    x="n_genes_by_counts",
    kind="reg",
)
plt.show()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)