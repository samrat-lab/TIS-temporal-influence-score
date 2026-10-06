function export_legacy_mats(scores, perm, cfg)
%EXPORT_LEGACY_MATS  Write the score matrices in the original cell-array layout.
%
%   EXPORT_LEGACY_MATS(scores, perm, cfg)
%
%   Earlier versions of this analysis stored results as {n_strain x 1} cell
%   arrays of string matrices whose first column is the TF name and whose
%   columns 2 ... n_tp+1 are the per-timepoint values:
%
%       score_mat_1st_term{k,1}  local influence L
%       score_mat_2nd_term{k,1}  network-propagated influence N
%       p_mat_1st_term{k,1}      permutation p-value of the combined score S
%
%   This function regenerates those variables so that downstream scripts
%   written against the original layout keep working. New code should use the
%   numeric matrices in `scores` / `perm` instead.
%
%   Note on p_mat: L and N are driven by the same random draw, so their
%   p-values are not independent; the combined p-value of S is written for
%   both terms, as in the original analysis.

score_mat_1st_term = cell(cfg.n_strain, 1);
score_mat_2nd_term = cell(cfg.n_strain, 1);
p_mat_1st_term     = cell(cfg.n_strain, 1);

for k = 1:cfg.n_strain
    score_mat_1st_term{k, 1} = [scores{k}.TF, string(scores{k}.L)];
    score_mat_2nd_term{k, 1} = [scores{k}.TF, string(scores{k}.N)];
    if ~isempty(perm)
        p = perm{k}.p_S;
        p(p == 0) = 1 / perm{k}.n_perm;
        p_mat_1st_term{k, 1} = [perm{k}.TF, string(p)];
    end
end

out_dir = fullfile(cfg.results_dir, 'legacy');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'score_mat_1st_term.mat'), 'score_mat_1st_term', '-v7.3');
save(fullfile(out_dir, 'score_mat_2nd_term.mat'), 'score_mat_2nd_term', '-v7.3');
if ~isempty(perm)
    save(fullfile(out_dir, 'p_mat_1st_term.mat'), 'p_mat_1st_term', '-v7.3');
end

fprintf('Legacy-format matrices written to %s\n', out_dir);
end
