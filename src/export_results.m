function export_results(scores, perm, cfg)
%EXPORT_RESULTS  Write score and p-value tables to cfg.results_dir as CSV.
%
%   EXPORT_RESULTS(scores, perm, cfg)
%
%   Files written, one set per strain:
%     results/scores/TIS_scores_<strain>.csv        TF x (L, N, S) per timepoint
%     results/permutation/TIS_pvalues_<strain>.csv  observed scores + p-values
%     results/permutation/TIS_significant_TFs_<strain>.csv
%     results/TIS_summary.csv                       one row per strain
%
%   p-values of exactly 0 are reported as 1/n_perm, the resolution limit of a
%   permutation test with n_perm draws.

score_dir = fullfile(cfg.results_dir, 'scores');
perm_dir  = fullfile(cfg.results_dir, 'permutation');
if ~exist(score_dir, 'dir'); mkdir(score_dir); end
if ~isempty(perm) && ~exist(perm_dir, 'dir'); mkdir(perm_dir); end

hrs = cfg.tp_hours;

for k = 1:cfg.n_strain
    strain = cfg.strain_labels{k};
    s      = scores{k};

    % ---- score matrix ---------------------------------------------------
    T = table(s.TF, 'VariableNames', {'TF'});
    for t = 1:cfg.n_tp
        T.(sprintf('L_%dh', hrs(t))) = s.L(:, t);
    end
    for t = 1:cfg.n_tp
        T.(sprintf('N_%dh', hrs(t))) = s.N(:, t);
    end
    for t = 1:cfg.n_tp
        T.(sprintf('S_%dh', hrs(t))) = s.S(:, t);
    end
    writetable(T, fullfile(score_dir, sprintf('TIS_scores_%s.csv', strain)));
end

if isempty(perm)
    fprintf('Exported score tables to %s\n', score_dir);
    return
end

summary = table('Size', [cfg.n_strain 5], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'strain', 'n_TF', 'n_sig_S', 'n_sig_L', 'n_sig_N'});

for k = 1:cfg.n_strain
    strain = cfg.strain_labels{k};
    r      = perm{k};

    p_S = r.p_S; p_S(p_S == 0) = 1 / r.n_perm;
    p_L = r.p_L; p_L(p_L == 0) = 1 / r.n_perm;
    p_N = r.p_N; p_N(p_N == 0) = 1 / r.n_perm;

    T = table(r.TF, any(r.sig_S, 2), any(r.sig_L, 2), any(r.sig_N, 2), ...
        'VariableNames', {'TF', 'sig_S_any', 'sig_L_any', 'sig_N_any'});
    for t = 1:cfg.n_tp
        h = hrs(t);
        T.(sprintf('S_%dh', h))      = r.S_obs(:, t);
        T.(sprintf('L_%dh', h))      = r.L_obs(:, t);
        T.(sprintf('N_%dh', h))      = r.N_obs(:, t);
        T.(sprintf('pS_%dh', h))     = p_S(:, t);
        T.(sprintf('pL_%dh', h))     = p_L(:, t);
        T.(sprintf('pN_%dh', h))     = p_N(:, t);
        T.(sprintf('Snull_%dh', h))  = r.S_null_mean(:, t);
    end
    writetable(T, fullfile(perm_dir, sprintf('TIS_pvalues_%s.csv', strain)));

    writetable(T(any(r.sig_S, 2), :), fullfile(perm_dir, ...
        sprintf('TIS_significant_TFs_%s.csv', strain)));

    summary.strain(k)  = string(strain);
    summary.n_TF(k)    = numel(r.TF);
    summary.n_sig_S(k) = sum(any(r.sig_S, 2));
    summary.n_sig_L(k) = sum(any(r.sig_L, 2));
    summary.n_sig_N(k) = sum(any(r.sig_N, 2));
end

writetable(summary, fullfile(cfg.results_dir, 'TIS_summary.csv'));

fprintf('\n=== Summary ===\n');
disp(summary);
fprintf('Exported tables to %s\n', cfg.results_dir);
end
