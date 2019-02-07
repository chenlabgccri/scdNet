function [p,adj_mat] = scdNet(data,GS_group,p_cutoff)

% The function analyzes a single-cell RNA-Seq dataset to identify gene
% pairs of which correlation coefficients changes significantly between
% sample groups. Within each group, scdNet calculates gene-gene interaction
% strengths. After adjusting the effects of sample sizes by Fisher
% transforamtion, scdNet tests the statistical significance of the change in
% gene-gene correlation between sample groups.
% 
% Inputs: DATA is a K-by-N numeric matrix (in double precision), which
% contains the log2 transformed expression values of K genes in N samples.
% NaN represent zero read count. Genes with multiple observations (e.g.,
% different probes and differnet splicing variants) should be merged.
% 
% GS_GROUP is a vector with length N yieled by GeneSetScoring. Values 1,
% -1, and 0 represent activation, inactivation, and no-call of a gene set,
% respectively, for a sample.
% 
% P_CUTOFF is the threshold on raw p-value to define significantly
% significant differential gene-gene interaction.
% 
% Outputs: P is a K-by-K symmetric matrix, denoting the significance of
% changes in gene-gene interaction across sample groups.
% 
% ADJ_MAT is a K-by-K symmetric adjacency matrix, of which ADJ_MAT(i,j)=1
% denotes a differential interaction pair of i and j (i.e., the i-j edge in
% the differential interaction network.
% 
% Reference: Yu-Chiao Chiu, Tzu-Hung Hsiao, Li-Ju Wang, Yidong Chen, and 
% Yu‐Hsuan Joni Shao. scdNet: a computational tool for single-cell
% differential network analysis. BMC systems biology. 2018;12(Suppl 8):124.

tic

num_gene = size(data,1);
groups = unique(GS_group);
num_group = length(groups);
corr = nan(num_gene,num_gene,num_group);
corr_fish = nan(size(corr));
num_sam_nonnan = nan(size(corr));
diff_corr_fish = nan(num_gene,num_gene);

for i=1:num_group
    tmp = data(:,GS_group==groups(i));
    [R P] = corrcoef(tmp','rows','pairwise');
    corr(:,:,i) = R;
    
    num_sam_nonnan(:,:,i) = double((~isnan(tmp)))*double((~isnan(tmp')));

    z1 = 0.5*log((1+corr(:,:,i))./(1-corr(:,:,i)));
    corr_fish(:,:,i) = sqrt(num_sam_nonnan(:,:,i)-3).*z1; 
end

for i=1:num_gene
    for j=1:num_gene
        diff_corr_fish(i,j) = max(abs(corr_fish(i,j,:))) - min(abs(corr_fish(i,j,:)));
    end
end

diff_corr_fish = -diff_corr_fish;

p = (0.5+erf(diff_corr_fish/2)-0.5*sign(diff_corr_fish).*erf(diff_corr_fish/2).*erf(diff_corr_fish/2));
adj_mat = double(p<p_cutoff);

time_used = toc;

disp(sprintf('\n\nSuccess! scdNet analysis was finished in %.2f seconds.\n',time_used));

