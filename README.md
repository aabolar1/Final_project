# Gene Expression Analysis Tool
This is a python program that performs multiple types of gene expression analysis on bulk mRNA sequencing data from  Wildtype and Resistnat Acute Myeloid Leukimic cells (specifcally Molm-13), which have been cultutred to be resistant to Venetoclax, a chemotherapeutic. 
## Analysis Types
- Differential Gene Expression Analysis (Diff exp)
- Principal Component Analysis (PCA)
- Gene Ontology Analysis (GO)
- Geneset enrichment Analysis (GSEA)

## Requirements: The following python packages are required to run this script
- pandas
- numpy
- seaborn
- matplotlib
- anndata
- rpy2
- scikit-learn
- goatools
- gseapy

## Usage
1. Download the Excel file containing the example geneset data.
2. Instantiate an RNASeqData object with the filename of the RNA-Seq data in Excel format.
3. Use the clean_data method to filter genes based on a minimum count threshold.
4. Instantiate a DiffExp object with the filtered gene data.
5. Use the run_deseq method to perform differential expression analysis on the filtered data.
6. Use the plot_heatmap method to plot a heatmap of the most differentially expressed genes.
7. Use the PCAplot method to plot the PCA for your filtered geneset.
8. Instantiate a GeneSetEnrichment object with the differential expression analysis dataframe
9. Use the gsea method to perfrom gene set enrichmen analysis for GO_Biological_Process_2021 & GO_Molecular_Function_2021
10. GO analysis unavailable at this time.
