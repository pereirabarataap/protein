import numpy as np
import pandas as pd
import seaborn as sb
from venn import venn
import matplotlib as mpl
from tqdm.notebook import tqdm
from Bio import ExPASy, SwissProt
from matplotlib import cm, pyplot as plt
from tqdm.notebook import tqdm as tqdm_notebook
from sklearn.metrics import davies_bouldin_score
from sklearn.cluster import AgglomerativeClustering as AC

class MplColorHelper:
    """
    Helper class to colour countplot
    """
    def __init__(self, cmap_name, start_val, stop_val):
        self.cmap_name = cmap_name
        self.cmap = plt.get_cmap(cmap_name)
        self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
    def get_rgb(self, val):
        rgb = self.scalarMap.to_rgba(val)[:-1]
        return rgb
    
def protein_correlation_clustering_go_count(protein_df, correlation_method="spearman"):
    """
    This function performs the following, in sequence:
        1. Compute correlation matrix of protein-protein from protein_df
        2. Hierarchical clustering of proteins using each correlation-row as feature vector 
            2.1 Ward-linkage is used to compute distances between proteins
            2.2 Optimal number of clusters is selected on lowest Davies Bouldin score
        3. Venn diagram of each cluster showing GO terms of respective cluster proteins
        4. Counts number of proteins in each cluster which contain cluster-specific GO terms
        
    protein_df ---> pandas.DataFrame
                    Column names must be protein accessions: str like "P05067"
                    Row values represent measures of respective protein: float
                    
        |   |     P98160 |   O00468-6 |   Q15149-4 |     P08238 |     P06733 |
        |:--|-----------:|-----------:|-----------:|-----------:|-----------:|
        | 0 | 1.5315e+09 | 2.5305e+09 | 9.382e+08  | 2.9324e+09 | 4.9473e+09 |
        | 1 | 1.8858e+09 | 3.1737e+09 | 8.2028e+08 | 3.4993e+09 | 3.6529e+09 |
        | 2 | 1.1927e+09 | 2.5913e+09 | 6.1717e+08 | 2.7219e+09 | 7.6858e+09 |
        | 3 | 1.4522e+09 | 2.6129e+09 | 9.3064e+08 | 3.2274e+09 | 3.6938e+09 |
        | 4 | 1.6498e+09 | 2.3794e+09 | 6.6648e+08 | 1.9489e+09 | 4.2946e+09 |
    
    correlation_method ---> str: {"pearson", "kendall", "spearman"}
                            default: "spearman"
                            
    """
    record_protein_gos = {}
    proteins = np.array(protein_df.columns.tolist())
    for protein in tqdm_notebook(proteins, desc="Retrieving protein records"):
        handle = ExPASy.get_sprot_raw(protein)
        record = SwissProt.read(handle)
        x_refs = record.cross_references
        record_protein_gos[protein] = set([x_ref for x_ref in x_refs if "GO" in x_ref])
        
    scores = []
    max_n_clusters = 6
    correlations = protein_df.corr(method=correlation_method).values
    tentative_clusters = np.linspace(2, max_n_clusters, max_n_clusters-1).astype(int)
    for n_clusters in tentative_clusters:
        ac = AC(n_clusters=n_clusters)
        ac.fit(correlations)
        score = davies_bouldin_score(correlations, ac.labels_)
        scores.append(score)

    best_n_clusters = tentative_clusters[np.array(scores).argmin()]
    ac = AC(n_clusters=best_n_clusters)
    ac.fit(correlations)

    cluster_gos = {}
    cluster_proteins = {}
    cluster_labels = np.unique(ac.labels_)

    protein_gos = {protein: {go[2] for go in record_protein_gos[protein]} for protein in proteins}
    # go[0]->"GO", [1]->"GO:0008150" [2]->description [3]->source
    for cluster in cluster_labels:
        cluster_gos[cluster] = set()
        cluster_proteins[cluster] = proteins[ac.labels_==cluster]
        for protein in cluster_proteins[cluster]:
            cluster_gos[cluster] = cluster_gos[cluster].union(protein_gos[protein])

    COL = MplColorHelper('viridis', 0, best_n_clusters-1)
    fig, axs = plt.subplots(1,2,dpi=100, figsize=(7*best_n_clusters,5))
    axs[0].plot(
        tentative_clusters,
        scores,
        color="C0",
        label="Davies Bouldin\n(lower is better)"
    )
    axs[0].set_xlabel("Number of clusters")
    axs[0].set_ylabel("Score value")
    axs[0].legend(title="Clustering Score")
    axs[0].set_title("Protein hierarchical clustering on protein correlations via ward-linkage")
    diagram = venn(cluster_gos, ax=axs[1], legend_loc=0, fontsize=9)
    axs[1].set_title("Venn diagram of unique GO terms per cluster")

    plt.tight_layout()
    plt.show()

    fig, axs = plt.subplots(1,best_n_clusters,dpi=100, figsize=(7*best_n_clusters,5), sharey=True)
    all_gos = set()
    for protein in proteins:
        all_gos = all_gos.union(protein_gos[protein])

    cluster_specific_gos = {}
    for cluster in cluster_labels:
        cluster_specific_gos[cluster] = set()
        for protein in cluster_proteins[cluster]:
            cluster_specific_gos[cluster] = cluster_specific_gos[cluster].union(protein_gos[protein])
        for other_cluster in cluster_labels:    
            if cluster != other_cluster:
                for protein in cluster_proteins[other_cluster]:
                    cluster_specific_gos[cluster] = cluster_specific_gos[cluster].difference(protein_gos[protein])

    #cluster_specific_protein_gos_df = pd.DataFrame(cluster_specific_protein_gos)
    cluster_specific_protein_gos = {}
    for cluster in cluster_labels:
        cluster_specific_protein_gos[cluster] = {}
        for protein in cluster_proteins[cluster]:
            cluster_specific_protein_gos[cluster][protein] = cluster_specific_gos[cluster].intersection(protein_gos[protein])

    cluster_specific_protein_gos_df = pd.DataFrame(cluster_specific_protein_gos)
    for cluster in cluster_labels:
        sb.countplot(
            x=cluster_specific_protein_gos_df[cluster].dropna().apply(lambda x: len(x)).sort_values(),
            ax=axs[cluster],
            color=COL.get_rgb(cluster),
            alpha=0.75
        )
        axs[cluster].set_ylabel("")
        axs[cluster].set_xlabel("Number of GOs")
        axs[cluster].set_title(
            "Cluster " + str(cluster) + " with " + 
            str(len(cluster_specific_gos[cluster])) + " cluster-specific GOs and" + 
            str(len(cluster_proteins[cluster])) + " proteins"
        )
    axs[0].set_ylabel("Number of proteins")
    plt.suptitle("How many proteins (y-axis) have this many GOs (x-axis) from the pool of cluster-specific GOs", fontsize=15)

    plt.tight_layout()
    plt.show()
