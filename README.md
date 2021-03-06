# scdNet
scdNet: a computational tool for single-cell differential network analysis


scdNet is devised to analyze differential gene regulatory networks associated with cellular states at the single cell level.
Briefly, scRNA-Seq data were preprocessed and normalized to eliminate inter-cell biases.

GeneSetScoring.m: For the per-cell state of a cellular function may not be fully represented by a gene or biomarker, we utilized a gene set analysis to quantify the activity of a function by summarizing the expressional abundance of all involved genes. Cells were
grouped according to function activities.

scdNet.m: Within each group, gene-gene correlation coefficients were calculated with an
exclusion of zero elements, and transformed to a sample size independent domain by the Fisher transformation to eliminate
sample size related biases. Normalized correlation coefficients were compared across groups of cells and the changes in correlation
were statistically tested in the Fisher domain. Significantly changed gene-gene pairs were merged into a differential
network.

Visualization and functional annotation analyses can be performed based on the outputs to realize the biological relevance of such dynamic network.

Authors:

Yu-Chiao Chiu1, Tzu-Hung Hsiao2, Li-Ju Wang1, Yidong Chen1,3, Yu‐Hsuan Joni Shao4

1Greehey Children's Cancer Research Institute, University of Texas Health Science Center at San Antonio, San Antonio, TX 78229, USA
2Department of Medical Research, Taichung Veterans General Hospital, Taichung 40705, Taiwan
3Department of Epidemiology and Biostatistics, University of Texas Health Science Center at San Antonio, San Antonio, TX 78229, USA
4Graduate Institute of Biomedical Informatics, College of Medical Science and Technology, Taipei Medical University, Taipei 10675, Taiwan

Manuscript of scdNet is accepted for oral presentation in The 2018 International Conference on Intelligent Biology and Medicine (ICIBM 2018) at Los Angeles, CA, USA, and published in BMC systems biology. 2018;12(Suppl 8):124.
