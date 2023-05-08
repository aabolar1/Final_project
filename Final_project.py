import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import rpy2.robjects as ro
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.conversion import localconverter
from sklearn.decomposition import PCA
deseq = importr('DESeq2')
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.base import download_go_basic_obo
import gseapy as gp

class RNASeqData:
    '''
    This class object represents the RNA seq data. It contains the name of the
    genes and the raw counts for each gene.
    '''

    def __init__(self, filename):
        self.filename = filename
        self.data = pd.read_excel(filename)
        self.gene_data_list = self.create_gene_data_list()

    def create_gene_data_list(self):
        '''
        This method creates a list of dictionaries, with each dictionary containing the gene information and raw counts
        for a single gene.
        
        Returns:
            A list of dictionaries, with each dictionary containing the gene information and raw counts for a single gene.
        '''
        gene_data_list = []
        for i, row in self.data.iterrows():
            counts = [row[f'Wildtype_{j}'] for j in range(1, 5)] + [row[f'Resistant_{j}'] for j in range(1, 5)]
            gene_data = {
                'gene_id': row['gene_id'],
                'gene_name': row['gene_name'],
                'counts': counts
            }
            gene_data_list.append(gene_data)
        return gene_data_list
 
    def clean_data(self, threshold):
        '''
        This method removes any genes with counts lower than the specified threshold.
        
        Returns:
            A list of dictionaries, with each dictionary containing the gene information and raw counts
            for a single gene.
        '''
        self.gene_data_list = [gene_data for gene_data in self.gene_data_list
                               if all(int(count) >= threshold for count in gene_data['counts'])]
        print(f"Total number of genes post filtering: {len(self.gene_data_list)}")
        return self.gene_data_list

class diff_Exp:
    '''
    This class object performs differential expression analysis on the filtered RNA-seq data.
    '''

    def __init__(self, gene_dataset):
        self.gene_data_list = gene_dataset
        self.data = self.create_data_frame()
        self.gene_name = self.get_gene_names()
        
    def get_gene_names(self):
        return [gene['gene_name'] for gene in self.gene_data_list]
    
    def create_data_frame(self):
        '''
        This method creates a Pandas DataFrame from the cleaned gene data list.

        Returns:
            A Pandas DataFrame containing gene information and raw counts for each sample.
        '''
        gene_ids = [gene['gene_id'] for gene in self.gene_data_list]
        gene_names = [gene['gene_name'] for gene in self.gene_data_list]
        counts = [gene['counts'] for gene in self.gene_data_list]
        
        data = pd.DataFrame(counts, columns=['Wildtype_1', 'Wildtype_2', 'Wildtype_3', 'Wildtype_4', 'Resistant_1', 'Resistant_2', 'Resistant_3', 'Resistant_4'])
        data.insert(0, 'gene_id', gene_ids)
        data.insert(1, 'gene_name', gene_names)
        data = data.infer_objects()
        return data
    
    def run_deseq(self, gene_names):
        '''
        Runs DESeq2 on the RNA-seq data and returns the results.
        '''
        # create sample_info DataFrame with correct number of conditions
        sample_names = self.data.columns[2:]
        conditions = [name.split('_')[0] for name in sample_names]
        sample_info = pd.DataFrame({'sample': sample_names, 
                                'condition': conditions})
        
        robjects.globalenv['sample_info'] = ro.conversion.py2rpy(sample_info)
        design_formula = robjects.Formula("~ condition")
        counts = self.data.iloc[:,2:].astype(int).values
        dds = deseq.DESeqDataSetFromMatrix(countData=counts,
                                           colData=sample_info,
                                           design=design_formula)
        dds = deseq.DESeq(dds)
        res = deseq.results(dds)
        res_r_df = robjects.r['as.data.frame'](res)
        with localconverter(ro.default_converter + pandas2ri.converter):
            res_pd = ro.conversion.rpy2py(res_r_df)
        res_pd.insert(0, 'gene_name', gene_names)
        return res_pd
    
    def plot_heatmap(self, res, top_percent=20):
        '''
        Plots a heatmap of the most differentially expressed genes.
        '''
        # Get the top x% most differentially expressed genes
        num_genes = int(len(res) * (top_percent / 100))
        top_genes = res.nlargest(num_genes, 'log2FoldChange')

        # Extract gene names and log2 fold change
        gene_names = top_genes['gene_name']
        log2fc = top_genes['log2FoldChange']

        # Create DataFrame for heatmap
        heatmap_data = pd.DataFrame({'Gene': gene_names, 'log2FoldChange': log2fc})
        heatmap_data = heatmap_data.set_index('Gene')
        heatmap_data = heatmap_data.sort_values('log2FoldChange')

        # Create color map
        cmap = plt.cm.RdBu_r

        # Create heatmap
        sns.set(font_scale=0.7)
        sns.heatmap(heatmap_data,xticklabels=False, yticklabels=True, cmap=cmap, linewidths=0.1, linecolor='gray')
        plt.xlabel('log2FoldChange')
        plt.ylabel('Gene')

        return plt
   
    def plot_volcano(self, res, top_percent=20):
        '''
        Plots a volcano plot of the differential expression results.
        '''
        # Set up plot
        plt.figure(figsize=(8, 8))
        plt.xlabel('log2(Fold Change)')
        plt.ylabel('-log10(p-value)')
        plt.title('Volcano Plot')

        # Plot points
        num_genes = int(len(res) * (top_percent / 100))
        top_genes = res.nlargest(num_genes, 'log2FoldChange')

        sig = top_genes[top_genes['padj'] <= 0.05]
        not_sig = top_genes[top_genes['padj'] > 0.05]
        plt.scatter(not_sig['log2FoldChange'], -np.log10(not_sig['padj']), color='grey', alpha=0.5,label='not significant')
        plt.scatter(sig['log2FoldChange'], -np.log10(sig['padj']), color='red', alpha=0.5, label='significant')

        # Add significance thresholds
        plt.axvline(x=-1, color='grey', linestyle='--')
        plt.axvline(x=1, color='grey', linestyle='--')
        plt.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
       
        # Add labels to significant genes
        for i, row in sig.iterrows():
            plt.text(row['log2FoldChange'], -np.log10(row['padj']), row['gene_name'], ha='center', va='bottom', fontsize=8)
        
        # Add legend
        plt.legend()
        
        return plt
    
class PCAPlot:
    '''
    A class for creating and visualizing PCA plots.
    '''
    def __init__(self, gene_dataset):
        self.gene_dataset = gene_dataset
        self.pca = PCA(n_components=2)

    def compute_transform(self):
        '''
        Fit a PCA model to the gene expression data and compute the
        transformed data in the PCA space.
        '''
        data = np.array([d['counts'] for d in self.gene_dataset])
        gene_names = [d['gene_name'] for d in self.gene_dataset]
        sample_labels = [name.split('_')[0] for name in gene_names]
        
        self.pca.fit(data)
        transformed_data = self.pca.transform(data)
        
        return transformed_data, sample_labels
    
    def plot(self):
        '''
        Create a scatter plot of the transformed data in the PCA space
        '''
        transformed_data, sample_labels = self.compute_transform()
        
        for label in set(sample_labels):
            index = [i for i, l in enumerate(sample_labels) if l == label]
            plt.scatter(transformed_data[index, 0], transformed_data[index, 1], alpha=0.5, label=label)
            
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCA Plot')
        plt.legend()
        plt.show()  
        
        return plt
    
class GeneSetEnrichment:
    '''
    A class that represents gene set enrichment analysis (GSEA).
    '''
    def __init__(self, filtered_gene_dataset):
        self.filtered_gene_dataset = filtered_gene_dataset
        
    def to_df(self):
        return self.data
        
    def run_gsea(self, gene_set, gene_sets_file, num_permutations=1000, figsize=(8, 8)):
        # perform GSEA using gseapy
        gs_res = gp.gsea(
            data=self.filtered_gene_dataset.to_df(), 
            gene_sets=self.gene_sets_file, 
            cls=gene_set,
            permutation_type='phenotype',
            permutation_num=num_permutations,
            method='signal_to_noise',
            verbose=False
        )
        
        # plot the GSEA results using matplotlib
        plt.figure(figsize=figsize)
        plt.plot(gs_res.es, gs_res.re, 'o', markersize=5)
        plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
        plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
        plt.xlabel('Enrichment score')
        plt.ylabel('Rank')
        plt.title(f'GSEA {gene_set} analysis')
        plt.show()         
    
        return plt
    
class GeneOntologyEnrichment:
    '''
    A class that performs Gene Ontology Enrichment Analysis.
    '''
    def __init__(self, filtered_data, species='human', pval_threshold=0.05, log_fc_threshold=1.0):
        self.filtered_data = filtered_data
        self.species = species
        self.pval_threshold = pval_threshold
        self.log_fc_threshold = log_fc_threshold
        self.geneids_pop = self.get_gene_ids(species)
        self.assoc = self.get_gene_associations(species)
        self.goea_results_all = self.run_go_analysis()
        
    def get_gene_ids(self, species):
        if species == 'human':
            geneids = read_ncbi_gene2go("gene2go", taxids=[9606])[0]
        return geneids
        
    def get_gene_associations(self, species):
        gene2go = read_ncbi_gene2go("gene2go", taxids=[9606])
        gene2go_dict = {}
        for i, j in gene2go.items():
            if j:
                gene2go_dict[i] = j
        return gene2go_dict
    
    def run_go_analysis(self):
        obo_dag = GODag("go-basic.obo")
        goeaobj = GOEnrichmentStudy(
            self.geneids_pop,
            self.assoc, 
            obo_dag, 
            propagate_counts=False,
            alpha=self.pval_threshold,
            log=None
        )
        geneids_study = self.filtered_gene_dataset.gene_names
        geneids_study = [gene for gene in geneids_study if gene in self.geneids_pop]
        return goeaobj.run_study(geneids_study)
    
    def plot_go_enrichment(self, num_enriched_terms=10):
        goea_results_sig = [r for r in self.goea_results_all if r.p_fdr_bh < self.pval_threshold and r.study_count > 2]
        goea_results_sig_sorted = sorted(goea_results_sig, key=lambda r: r.p_fdr_bh)
        enriched_terms = [r.GO for r in goea_results_sig_sorted][:num_enriched_terms]
        
if __name__ == '__main__':
    # Load and clean RNA sequencing data
    rna_seq_data = RNASeqData('gene_fpkm.xlsm')
    filtered_data = rna_seq_data.clean_data(1)
    
    # Run differential expression analysis
    diff = diff_Exp(filtered_data)
    res_pd = diff.run_deseq(diff.get_gene_names())
    
    # Plot heatmap and volcano plot
    heatmap = diff.plot_heatmap(res_pd)
    heatmap.savefig('heatmap.png')
    plt.show()
    
    volcano = diff.plot_volcano(res_pd)
    volcano.savefig('volcano.png')
    plt.show()
    
    # Create PCA plot of gene expression data
    pca_plot = PCAPlot(filtered_data)
    pca = pca_plot.plot()
    pca.savefig('pca.png')
    plt.show()
    
    # Perform GeneSetEnrichment
    #enrichment = gsea.run_gsea('h.all.v2023.1.Hs.symbols.gmt.txt', '/Users/alliabolarin/Documents/Python/h.all.v2023.1.Hs.symbols.gmt.txt')
    
     
    # Perform GO analysis
    #go_analysis = GoAnalysis(filtered_data, res_pd)
    #go_analysis.run_go_analysis()
    #go_analysis.plot_go_enrichment()
