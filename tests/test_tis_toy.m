function test_tis_toy()
%TEST_TIS_TOY  Self-contained unit test of the TIS implementation.
%
%   Run from the repository root:   test_tis_toy
%
%   A four-timepoint toy network with hand-computed scores exercises every
%   branch of the recursion: the terminal transition (N = 0), a transition
%   with both a non-TF target and a downstream TF, and the layered network /
%   reach / score-matrix assembly.
%
%   Toy design
%     t1  TFs {A}      proteins {A, P1}
%     t2  TFs {B}      proteins {B, P2, P3}
%     t3  TFs {C}      proteins {C, P4, P5}
%     t4  TFs {}       proteins {P6, P7}
%     TFLink edges  A->P2, A->B, A->P9, B->P4, B->C, C->P6, C->P7
%                   (A->P9 is never measured and must be ignored)
%
%   Hand calculation
%     i = 3 (terminal) : P_4 = {P6,P7} so L_C = 2/2 = 1,   N_C = 0,   S_C = 1.0
%     i = 2            : P_3 = {P4,P5} so L_B = 1/2 = 0.5
%                        N_B = S_C / S_C = 1.0            S_B = 1.5
%     i = 1            : P_2 = {P2,P3} so L_A = 1/2 = 0.5
%                        N_A = S_B / S_B = 1.0            S_A = 1.5

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'config'));

cfg = tis_config('strain_labels', {'TOY'}, 'tp_hours', [1 2 3 4], ...
                 'verbose', false, 'n_perm', 50, 'rng_seed', 1);

tol = 1e-12;

% ---- synthetic inputs ---------------------------------------------------
prot = {["A";"P1"], ["B";"P2";"P3"], ["C";"P4";"P5"], ["P6";"P7"]};
tfs  = {"A",        "B",             "C",             strings(0,1)};

data_info = cell(cfg.n_tp, 2);
for t = 1:cfg.n_tp
    blk = cell(1, 6);
    blk{1, cfg.col_master}      = prot{t};
    blk{1, cfg.col_newly_synth} = prot{t};
    blk{1, cfg.col_up}          = prot{t};
    blk{1, cfg.col_down}        = strings(0, 1);
    blk{1, cfg.col_densp}       = prot{t};
    blk{1, cfg.col_tf}          = string(tfs{t}(:));
    data_info{t, 1} = cfg.tp_hours(t);
    data_info{t, 2} = blk;
end

tflink = ["A" "P2"; "A" "B"; "A" "P9"; "B" "P4"; "B" "C"; "C" "P6"; "C" "P7"];

% ---- pipeline ------------------------------------------------------------
TF_target_info = build_tf_target_info(data_info, tflink, cfg);
nodes          = build_node_sets(data_info, cfg);
tis            = compute_tis_scores(nodes, TF_target_info, cfg);
net            = build_temporal_network(nodes, TF_target_info, cfg);
scores         = assemble_score_matrices(tis, TF_target_info, cfg);

% ---- 1. TF universe and retained TFLink edges ---------------------------
assert(isequal(TF_target_info{1,1}, ["A";"B";"C"]), 'TF universe is wrong.');
assert(size(TF_target_info{1,2}, 1) == 7, 'All 7 TFLink edges should be retained.');

% ---- 2. recursive scores -------------------------------------------------
check(tis{1}.L{3}, 1.0, tol, 'L at terminal transition (C)');
check(tis{1}.N{3}, 0.0, tol, 'N at terminal transition (C)');
check(tis{1}.L{2}, 0.5, tol, 'L at i=2 (B)');
check(tis{1}.N{2}, 1.0, tol, 'N at i=2 (B)');
check(tis{1}.S{2}, 1.5, tol, 'S at i=2 (B)');
check(tis{1}.L{1}, 0.5, tol, 'L at i=1 (A)');
check(tis{1}.N{1}, 1.0, tol, 'N at i=1 (A)');
check(tis{1}.S{1}, 1.5, tol, 'S at i=1 (A)');
assert(isempty(tis{1}.S{4}) || all(tis{1}.S{4} == 0), ...
    'No transition starts at the last timepoint.');

% ---- 3. unmeasured target A->P9 must not contribute ---------------------
assert(abs(tis{1}.L{1} - 0.5) < tol, 'Unmeasured target leaked into L.');

% ---- 4. temporal network -------------------------------------------------
expected_edges = sortrows(["A_1" "B_2"; "A_1" "P2_2"; "B_2" "C_3"; ...
                           "B_2" "P4_3"; "C_3" "P6_4"; "C_3" "P7_4"]);
assert(isequal(sortrows(net{1}.edges), expected_edges), 'Temporal edge set is wrong.');

% ---- 5. score matrices ---------------------------------------------------
expL = [0.5 0 0 0; 0 0.5 0 0; 0 0 1.0 0];
expN = [1.0 0 0 0; 0 1.0 0 0; 0 0 0   0];
assert(isequal(scores{1}.TF, ["A";"B";"C"]), 'Score matrix row labels are wrong.');
assert(max(abs(scores{1}.L(:) - expL(:))) < tol, 'L matrix is wrong.');
assert(max(abs(scores{1}.N(:) - expN(:))) < tol, 'N matrix is wrong.');
assert(max(abs(scores{1}.S(:) - (expL(:) + expN(:)))) < tol, 'S matrix is wrong.');
assert(all(scores{1}.S(:, cfg.n_tp) == 0), ...
    'No score is defined at the final timepoint.');

% ---- 6. permutation test runs and returns valid p-values ----------------
perm = run_permutation_test(nodes, TF_target_info, scores, cfg);
p    = perm{1}.p_S;
assert(isequal(size(p), [3 cfg.n_tp]), 'p-value matrix has the wrong size.');
assert(all(p(:) >= 0 & p(:) <= 1), 'p-values must lie in [0,1].');
assert(all(p(:, cfg.n_tp) == 1), 'No transition starts at the last timepoint.');

fprintf('test_tis_toy: all checks passed.\n');
end

% -------------------------------------------------------------------------
function check(actual, expected, tol, label)
assert(isscalar(actual), '%s: expected a scalar, got %d values.', label, numel(actual));
assert(abs(actual - expected) < tol, ...
    '%s: expected %.6f, got %.6f.', label, expected, actual);
end
