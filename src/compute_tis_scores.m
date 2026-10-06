function tis = compute_tis_scores(nodes, TF_target_info, cfg)
%COMPUTE_TIS_SCORES  Recursive temporal influence score (Eqs. 8-12).
%
%   tis = COMPUTE_TIS_SCORES(nodes, TF_target_info, cfg)
%
%   For every strain k and every transition i -> i+1 (i = 1 ... n_tp-1) the
%   score of a transcription factor t detected at timepoint i is
%
%       S^(i)_t = L^(i->i+1)_t + N^(i->i+1)_t                            (Eq. 8)
%
%   with the LOCAL INFLUENCE (1st term, Eq. 10)
%
%       L^(i->i+1)_t = | N_target(t) | / | P_{i+1} |
%
%   where N_target(t) are the TFLink targets of t that are measured at
%   timepoint i+1 and are NOT themselves transcription factors, and P_{i+1}
%   is the complete set of newly synthesized non-TF proteins at i+1;
%
%   and the NETWORK-PROPAGATED INFLUENCE (2nd term, Eq. 11)
%
%       N^(i->i+1)_t = sum_{u in N_TF(t)} S^(i+1)_u  /  sum_{v in T_{i+1}} S^(i+1)_v
%
%   where N_TF(t) are the TFLink targets of t that are measured at i+1 and
%   are themselves transcription factors.
%
%   The recursion is evaluated BACKWARDS in time. At the terminal transition
%   (i = n_tp-1) no downstream TF layer with defined scores exists, so
%   N = 0 and S = L (Eq. 9).
%
%   Returned structure, per strain k:
%     tis{k}.TF{i}  [n_i x 1] string  - TFs detected at timepoint i
%     tis{k}.L{i}   [n_i x 1] double  - 1st term
%     tis{k}.N{i}   [n_i x 1] double  - 2nd term
%     tis{k}.S{i}   [n_i x 1] double  - L + N
%   Entries for i = n_tp are empty: no transition starts at the last timepoint.

tis = cell(cfg.n_strain, 1);
nT  = cfg.n_tp;

for k = 1:cfg.n_strain
    src = TF_target_info{k, 2}(:, 1);
    dst = TF_target_info{k, 2}(:, 2);

    s    = struct();
    s.TF = cell(nT, 1);
    s.L  = cell(nT, 1);
    s.N  = cell(nT, 1);
    s.S  = cell(nT, 1);
    for i = 1:nT
        s.TF{i} = nodes{k}.tf{i};
        s.L{i}  = zeros(numel(s.TF{i}), 1);
        s.N{i}  = zeros(numel(s.TF{i}), 1);
        s.S{i}  = zeros(numel(s.TF{i}), 1);
    end

    % Backward recursion over transitions i -> i+1.
    for i = (nT - 1):-1:1
        tf_i      = s.TF{i};
        prot_next = nodes{k}.prot{i+1};
        tf_next   = nodes{k}.tf{i+1};
        nonTF_next= nodes{k}.nonTF{i+1};

        n_total_target = numel(nonTF_next);          % | P_{i+1} |
        is_terminal    = (i == nT - 1);

        if ~is_terminal
            S_next    = s.S{i+1};
            denom_TF  = sum(S_next);                 % sum over all TFs at i+1
        end

        L = zeros(numel(tf_i), 1);
        N = zeros(numel(tf_i), 1);

        for j = 1:numel(tf_i)
            % ---- forward temporal neighbourhood  N_t  (Eq. 5) ----------
            tgt = unique(dst(src == tf_i(j)));
            tgt = intersect(tgt, prot_next);         % only measured targets

            % ---- 1st term: non-TF targets (Eqs. 6, 10) ------------------
            if n_total_target > 0
                L(j) = numel(setdiff(tgt, tf_next)) / n_total_target;
            end

            % ---- 2nd term: downstream TF targets (Eqs. 7, 11) -----------
            if ~is_terminal
                [~, ~, idx_next] = intersect(tgt, tf_next);
                if ~isempty(idx_next) && denom_TF > 0
                    N(j) = sum(S_next(idx_next)) / denom_TF;
                end
            end
        end

        s.L{i} = L;
        s.N{i} = N;
        s.S{i} = L + N;
    end

    tis{k} = s;

    if cfg.verbose
        nz = sum(cellfun(@(x) sum(x > 0), s.S));
        fprintf('%-8s : %d TF-timepoint pairs with S > 0\n', ...
            cfg.strain_labels{k}, nz);
    end
end
if cfg.verbose; fprintf('\n'); end
end
