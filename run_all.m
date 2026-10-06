function out = run_all(varargin)
%RUN_ALL  End-to-end driver for the Temporal Influence Score (TIS) pipeline.
%
%   out = RUN_ALL()
%   out = RUN_ALL('n_perm', 200)
%   out = RUN_ALL('run_permutation', false)
%
%   Reproduces, for every strain and every temporal transition:
%     * the directed temporal TF -> target network (Eqs. 1-7 of the paper),
%     * the local influence term        L  (1st term, Eq. 10),
%     * the network-propagated term     N  (2nd term, Eq. 11),
%     * the temporal influence score    S = L + N (Eq. 12),
%     * permutation p-values for S, L and N (Eq. 13).
%
%   All name-value options are forwarded to TIS_CONFIG (see config/tis_config.m).
%
%   Inputs expected in ./data :
%     data_info.mat        - proteomics data structure (see data/README.md)
%     Final_TF_Target.mat  - TFLink TF -> target edge list (source, target)
%
%   Outputs written to ./results .

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
addpath(fullfile(here, 'config'));

cfg = tis_config(varargin{:});
if ~exist(cfg.results_dir, 'dir'); mkdir(cfg.results_dir); end

fprintf('===============================================================\n');
fprintf(' Temporal Influence Score (TIS) pipeline\n');
fprintf('===============================================================\n');
fprintf(' Strains    : %s\n', strjoin(cfg.strain_labels, ', '));
fprintf(' Timepoints : %s h\n', strjoin(string(cfg.tp_hours), ', '));
fprintf(' Results    : %s\n\n', cfg.results_dir);

t0 = tic;

% ---- 1. inputs ---------------------------------------------------------
[data_info, tflink] = load_inputs(cfg);

% ---- 2. strain-specific TF universe and TFLink sub-network -------------
TF_target_info = build_tf_target_info(data_info, tflink, cfg);

% ---- 3. per-timepoint protein / TF node sets ---------------------------
nodes = build_node_sets(data_info, cfg);

% ---- 4. recursive TIS (1st term L, 2nd term N, total S) ----------------
tis = compute_tis_scores(nodes, TF_target_info, cfg);

% ---- 5. directed temporal network (Section 2.3) ------------------------
% Returned for inspection and export; the score itself realises the same
% edges inside compute_tis_scores.
net = build_temporal_network(nodes, TF_target_info, cfg);

% ---- 6. TF x timepoint score matrices ----------------------------------
scores = assemble_score_matrices(tis, TF_target_info, cfg);

% ---- 7. permutation test ------------------------------------------------
if cfg.run_permutation
    perm = run_permutation_test(nodes, TF_target_info, scores, cfg);
else
    perm = [];
    fprintf('Permutation test skipped (run_permutation = false).\n\n');
end

% ---- 8. export ----------------------------------------------------------
export_results(scores, perm, cfg);
export_legacy_mats(scores, perm, cfg);

out               = struct();
out.cfg           = cfg;
out.TF_target_info= TF_target_info;
out.nodes         = nodes;
out.tis           = tis;
out.net           = net;
out.scores        = scores;
out.perm          = perm;

save(fullfile(cfg.results_dir, 'tis_pipeline.mat'), '-struct', 'out', '-v7.3');
fprintf('Saved workspace: %s\n', fullfile(cfg.results_dir, 'tis_pipeline.mat'));
fprintf('Total runtime: %.1f s\n', toc(t0));
end
