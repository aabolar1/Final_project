# Final project for Software Carpentry Spring 2023
This is a python program that performs multiple types of gene expression analysis on bulk mRNA sequencing data. The samples Wildtype and Resistnat Acute Myeloid Leukimic cells, specifcally Molm-13, which have been cultutred to be resistant to Venetoclax, a chemotherapeutic.
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
- Instantiate an RNASeqData object with the filename of the RNA-Seq data in Excel format
- Use the clean_data method to filter genes based on a minimum count threshold
- Instantiate a diff_Exp object with the filtered gene data
- Use the run_deseq method to perform differential expression analysis on the filtered data
- Use the plot_heatmap method to plot a heatmap of the most differentially expressed genes.
