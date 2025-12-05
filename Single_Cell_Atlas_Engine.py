import scanpy as sc
import matplotlib.pyplot as plt
import os

# Suppress tensorflow warnings if present
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

class ProjectX:
    def pipeline(self, adata):
        # 1. Load & QC
        adata = self.load_and_clean_data(adata)
        if adata is None: return

        # 2. Normalize & Select Features
        adata = self.normalizer_and_feature_selector(adata)

        # 3. Dimensionality Reduction (PCA -> UMAP)
        adata = self.dimensionality_reduction(adata)

        # 4. Clustering & Marker Identification
        adata = self.clustering_and_annotation(adata)
        
        return adata
    
    def load_and_clean_data(self, file_path : str):
        try:
            # Load the data using scanpy
            file = sc.datasets.pbmc3k()
            file.var["mito"] = file.var_names.str.startswith("MT-")

            print(f"Original data shape: {file.shape}") # Shape of the data matrix (2700 x 32738)

            # Calculate QC metrics
            sc.pp.calculate_qc_metrics(file, qc_vars = ["mito"], inplace=True)

            # Visualize the seaborn plot
            # sns.jointplot(data=file.obs, x = "log1p_total_counts_mito", y="log1p_total_counts", kind="hex")
            # plt.show()

            print("Filtering under progress...")

            # Filter the cells
            # Remove cells having mitochondrial reads > 5% (Dead cells)
            filtered_file = file[file.obs["total_counts_mito"] / file.obs["total_counts"] < 0.05].copy()

            # Remove cells with < 200 genes
            sc.pp.filter_cells(filtered_file, min_genes = 200, inplace = True)

            # Remove cells with > 2500 genes
            sc.pp.filter_cells(filtered_file,  max_genes = 2500, inplace = True)

            print(f"Filtered data shape: {filtered_file.shape}") # Shape of the filtered file (2638 x 32738)
            return filtered_file

        except FileNotFoundError:
            print(f"File not found at {file_path}")
            return None
        
    def normalizer_and_feature_selector(self, filtered_file):
        try:
            # Normalize the data (CPM)
            sc.pp.normalize_total(filtered_file, target_sum = 1e4, inplace = True)
            sc.pp.log1p(filtered_file)

            # Filter the matrix to keep only highly variable genes
            sc.pp.highly_variable_genes(filtered_file, min_mean=0.0125, max_mean=3, min_disp=0.5, inplace = True)

            # Plot the variable genes
            sc.pl.highly_variable_genes(filtered_file, show=False)
            plt.title("Highly Variable Genes")
            plt.show()

            filtered_file.raw = filtered_file

            filtered_file = filtered_file[:, filtered_file.var.highly_variable]

            filtered_file = filtered_file.copy()

            sc.pp.scale(filtered_file, max_value = 10, zero_center=False)

            print(f"Shape after feature selection: {filtered_file.shape}") # 2638 x 2013
            return filtered_file

        except FileNotFoundError:
            print(f"Error: Filtered file not found at {filtered_file}")
            return None
        
    def dimensionality_reduction(self, filtered_file):
        try:
            # Run PCA
            sc.pp.pca(filtered_file, svd_solver = "arpack")

            # Plot Variance Ratio (Elbow Plot) to see how many PCs matter
            sc.pl.pca_variance_ratio(filtered_file, log=True, show=False)
            plt.show()

            # Run UMAP after PCA to reduce dimensionality to 2
            sc.pp.neighbors(filtered_file, n_neighbors=10, n_pcs=40)
            sc.tl.umap(filtered_file)

            # Visualize UMAP embeddings
            sc.pl.umap(filtered_file, color=['CST3'], show=False)
            plt.title("UMAP: CST3 Expression")
            plt.show()
            return filtered_file

        except FileNotFoundError:
            print(f"Error: Filtered file not found at {filtered_file}")
            return None
        
    def clustering_and_annotation(self, filtered_file):
        try:
            sc.tl.leiden(filtered_file, resolution = 0.5, flavor = "igraph", n_iterations = 2, directed = False)
            sc.tl.rank_genes_groups(filtered_file, 'leiden', method='t-test')

            # Plot the ranked genes (Top 5 genes for the first few clusters)
            sc.pl.rank_genes_groups(filtered_file, n_genes=5, sharey=False, show=False)
            plt.show()

            # Plot the leiden graph
            sc.pl.umap(filtered_file, color=['leiden'], legend_loc='on data', title="Leiden Clusters", show=False)
            plt.show()
            return filtered_file

        except FileNotFoundError:
            print(f"Error: File not found at {filtered_file}")
            return None
        
sol = ProjectX()

# Following data file should be downloaded from the official website of 10x Genomics (https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
data_file : str = 'E:\\Project X\\pbmc3k_filtered_gene_bc_matrices\\filtered_gene_bc_matrices\\hg19'

sol.pipeline(data_file)

# ------------- Project Completed! ---------------
