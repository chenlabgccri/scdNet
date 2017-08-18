# scDNA
Single-cell RNA sequencing-based differential network analysis


scDNA is devised to analyze differential gene regulatory networks associated with cellular states at the single cell level.
Briefly, scRNA-Seq data were preprocessed and normalized to eliminate inter-cell biases. For the per-cell state of a
cellular function may not be fully represented by a gene or biomarker, we utilized a gene set analysis to quantify the
activity of a function by summarizing the expressional abundance of all involved genes. Cells were grouped according to
function activities. Within each group, gene-gene correlation coefficients were calculated with an exclusion of zeros
elements, and transformed to a sample size independent domain by the Fisher transformation to eliminate sample size
related biases. Normalized correlation coefficients were compared across groups of cells and the changes in correlation
were statistically tested in the Fisher domain. Significantly changed gene-gene pairs were merged into a differential
network. Visualization and functional annotation analyses were performed to realize the biological relevance of such
dynamic network.

Authors:

Yu-Chiao Chiu1, Yidong Chen1,2

1Greehey Children’s Cancer Research Institute
2Department of Epidemiology and Biostatistics
University of Texas Health Science Center at San Antonio
San Antonio, TX 78229, USA

Li-Ju Wang3, Tzu-Hung Hsiao3

3Department of Medical Research
Taichung Veterans General Hospital
Taichung 40705, Taiwan

Manuscript of scDNA is under consideration for a presentation in The IEEE International Conference on Bioinformatics and
Biomedicine (BIBM) 2017 at Kansas City, MO, USA.
