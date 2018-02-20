function [GS_score GS_group] = GeneSetScoring(data,data_gene,genes_OI,flag_Ztrans)

% The function performs a gene set analysis of a list of genes of interest
% in a gene expression profile. GeneSetScoring calculates the per-sample
% score of a gene set and statistically tests whether a score is
% significantly high (activation of the functions) or low (inactivation)
% and returns a vector denoting the groups of samples with respect to the
% function.
% 
% Inputs: DATA is a K-by-N numeric matrix (in double precision), which
% contains the expression profiles of K genes in N samples. DATA should not
% contain NaNs. Genes with multiple observations (e.g., different probes
% and differnet splicing variants) should be merged.
% 
% DATA_GENE is a cell vector with length K, containing gene symbols of
% genes in DATA. DATA_GENE should not contain numeric or empty elements.
% 
% GENES_OI is a cell vector of gene symbols of interest (i.e., members of a
% gene set).
% 
% FLAG_ZTRANS is set as 1 to perform z-transformation to each gene of DATA.
% If input DATA is already z-transformed, set FLAG_ZTRANS to 0.
% 
% Outputs: GS_SCORE is a vector with length N, containing per-sample gene
% set scores.
% 
% GS_GROUP is a vector with length N. Values 1, -1, and 0
% represent activation, inactivation, and no-call of a gene set,
% respectively, for a sample, with a alpha cutoff at 0.05.
% 
% Reference: Yu-Chiao Chiu, Tzu-Hung Hsiao, Li-Ju Wang, Yidong Chen, and 
% Yu‐Hsuan Joni Shao. scDNA: a fast and comprehensive tool for single-cell
% differential network analysis. Manuscript is under consideration for a
% presentation in The 2018 International Conference on Intelligent Biology
% and Medicine (ICIBM 2018) at Los Angeles, California, USA.

if flag_Ztrans==1
    for i=1:size(data,1)
        data(i,:) = (data(i,:)-mean(data(i,:)))/std(data(i,:));
    end
end

[tf loc] = ismember(data_gene,genes_OI);

GS_score = mean(data(tf~=0),:);
GS_p = 2*(1-normcdf(abs(GS_score),0,1/sqrt(sum(tf))));

GS_group = zeros(1,length(GS_score));
GS_group(GS_p<0.05) = sign(GS_score(GS_p<0.05));

