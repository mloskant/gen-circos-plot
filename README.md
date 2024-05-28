# gen-circos-plot
The program gen-circos-plot is used to visualize 2D miRNA folding and protein binding sites. It uses miRNA coordinates (BED file), the mapped binding sites (BED file) and the genome (FASTA file) to create a CIRCOS plot of the folded miRNA with binding sites marked. The package pyCircos is used for this.

# Usage
## Installation of pyCircos
pip install python-circos
## Command line flags are:
- -g or --genes : The genes in BED format
- -b or --bindings : The binding sites of proteins in BED format
- -fi or --fasta : The chromosomal sequences as a fasta file
- -o or --output : Name of the folder in which all pngs are saved

## Example: 
python3 gencircosplot.py --genes sample_data/example_data_genes.bed --bindings sample_data/example_data_binding_sites.bed --fasta sample_data/example_data_TAIR10_chr_all.fas -o circoplots
