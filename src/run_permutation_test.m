function perm = run_permutation_test(nodes, TF_target_info, scores, cfg)
%RUN_PERMUTATION_TEST  Empirical significance of the temporal influence score.
%
%   perm = RUN_PERMUTATION_TEST(nodes, TF_target_info, scores, cfg)
%
%   Null model (Eq. 13). For each strain and each transition i -> i+1 the set
%   of newly synthesized proteins observed at timepoint i+1 is replaced by an
%   equally sized sample drawn WITHOUT replacement from the union of proteins
%   detected across all timepoints in that strain. The TF -> target topology
%   is left untouched. L and N are recomputed on the permuted node set, giving
%   a null score S_null, and the empirical p-value is
%
%       p^(i)_t = (1/n) * sum_k 1[ S^(i)_{t,null,k} >= S^(i)_t ]
%
%   with n = cfg.n_perm. Separate p-values are also accumulated for the L and
%   N components. A TF is called significant if p < cfg.alpha at >= 1
%   transition.
%
%   Two properties of the null model are worth stating explicitly:
%     1. TF / non-TF status of a sampled protein is taken from the strain-wide
%        TF universe (TF_target_info{k,1}), not from the per-timepoint TF list,
%        because a randomly drawn protein has no timepoint of its own.
%     2. The propagated term uses the OBSERVED downstream scores S^(i+1) as
%        edge weights; only the membership of the downstream layer is
%        randomised. The test therefore asks whether a TF reaches an unusually
%        influential part of the next layer, not whether the whole recursive
%        cascade could arise by chance.
%
%   Returned structure, per strain k:
%     perm{k}.TF, .L_obs, .N_obs, .S_obs
%     perm{k}.p_S, .p_L, .p_N      [nTF x n_tp] empirical p-values
%     perm{k}.sig_S, .sig_L, .sig_N  logical, p < alpha
%     perm{k}.S_null_mean          [nTF x n_tp] mean null score
%     perm{k}.sig_TFs              TFs significant at >= 1 transition
%     perm{k}.n_perm, .alpha

perm = cell(cfg.n_strain, 1);
nT   = cfg.n_tp;

rng(cfg.rng_seed, 'twister');

fprintf('=== Permutation test (%d permutations, alpha = %.2f) ===\n', ...
    cfg.n_perm, cfg.alpha);

for k = 1:cfg.n_strain
    t_start = tic;
    strain  = cfg.strain_labels{k};

    % ---- protein pool: union of all proteins detected in this strain ----
    pool = strings(0, 1);
    for t = 1:nT
        pool = union(pool, nodes{k}.prot{t});
    end
    pool(pool == "") = [];
    n_pool = numel(pool);

    TF_list = TF_target_info{k, 1};
    nTF     = numel(TF_list);

    % ---- which pool entries are transcription factors? -------------------
    [tf_in_pool, tf_pool_idx] = ismember(TF_list, pool);
    is_TF_in_pool              = false(n_pool, 1);
    is_TF_in_pool(tf_pool_idx(tf_in_pool)) = true;

    % ---- TF -> target incidence matrices over the pool -------------------
    src = TF_target_info{k, 2}(:, 1);
    dst = TF_target_info{k, 2}(:, 2);
    [dst_in, dst_idx] = ismember(dst, pool);

    rows_non = []; cols_non = [];
    rows_tf  = []; cols_tf  = [];
    for ti = 1:nTF
        sel = (src == TF_list(ti)) & dst_in;
        di  = unique(dst_idx(sel));
        if isempty(di); continue; end
        dn = di(~is_TF_in_pool(di));
        dt = di( is_TF_in_pool(di));
        rows_non = [rows_non; repmat(ti, numel(dn), 1)];   %#ok<AGROW>
        cols_non = [cols_non; dn];                         %#ok<AGROW>
        rows_tf  = [rows_tf;  repmat(ti, numel(dt), 1)];   %#ok<AGROW>
        cols_tf  = [cols_tf;  dt];                         %#ok<AGROW>
    end
    M_non = sparse(rows_non, cols_non, 1, nTF, n_pool);   % non-TF targets
    M_tf  = sparse(rows_tf,  cols_tf,  1, nTF, n_pool);   % downstream TF targets

    % ---- observed quantities ---------------------------------------------
    L_obs = scores{k}.L;
    N_obs = scores{k}.N;
    S_obs = scores{k}.S;

    TF_present = false(nTF, nT);
    for t = 1:nT
        TF_present(:, t) = ismember(TF_list, nodes{k}.tf{t});
    end

    % Pool-indexed vector of downstream scores, one column per timepoint.
    W = zeros(n_pool, nT);
    valid = tf_in_pool;
    W(tf_pool_idx(valid), :) = S_obs(valid, :);

    n_obs = zeros(nT, 1);
    for t = 1:nT
        n_obs(t) = numel(nodes{k}.prot{t});
    end

    fprintf('  %-8s pool = %d proteins | TFs = %d\n', strain, n_pool, nTF);

    % ---- pre-generate all random draws -----------------------------------
    draws = cell(nT - 1, 1);
    for t = 1:(nT - 1)
        n_draw = min(n_obs(t + 1), n_pool);
        d      = zeros(n_draw, cfg.n_perm, 'uint32');
        for p = 1:cfg.n_perm
            d(:, p) = randperm(n_pool, n_draw);
        end
        draws{t} = d;
    end

    % ---- permutation loop --------------------------------------------------
    cnt_S = zeros(nTF, nT);
    cnt_L = zeros(nTF, nT);
    cnt_N = zeros(nTF, nT);
    sum_S = zeros(nTF, nT);

    denom_TF = zeros(nT, 1);
    for t = 1:nT
        col         = S_obs(:, t);
        denom_TF(t) = sum(col(col > 0));
    end

    for p = 1:cfg.n_perm
        for t = 1:(nT - 1)
            idx          = double(draws{t}(:, p));
            in_perm      = false(n_pool, 1);
            in_perm(idx) = true;

            Ni_p = sum(in_perm & ~is_TF_in_pool);
            if Ni_p == 0; continue; end

            L_p = (M_non * double(in_perm)) / Ni_p;

            if t < (nT - 1) && denom_TF(t + 1) > 0
                N_p = (M_tf * (double(in_perm) .* W(:, t + 1))) / denom_TF(t + 1);
            else
                N_p = zeros(nTF, 1);
            end

            S_p  = L_p + N_p;
            mask = TF_present(:, t);

            sum_S(mask, t) = sum_S(mask, t) + S_p(mask);
            cnt_S(mask, t) = cnt_S(mask, t) + (S_p(mask) >= S_obs(mask, t));
            cnt_L(mask, t) = cnt_L(mask, t) + (L_p(mask) >= L_obs(mask, t));
            cnt_N(mask, t) = cnt_N(mask, t) + (N_p(mask) >= N_obs(mask, t));
        end
        if cfg.verbose && mod(p, 100) == 0
            fprintf('    %4d / %d permutations (%.1f s)\n', p, cfg.n_perm, toc(t_start));
        end
    end

    % ---- p-values ----------------------------------------------------------
    p_S = cnt_S / cfg.n_perm;
    p_L = cnt_L / cfg.n_perm;
    p_N = cnt_N / cfg.n_perm;

    % A score that is absent or exactly zero cannot be called significant.
    p_S(~TF_present) = 1; p_L(~TF_present) = 1; p_N(~TF_present) = 1;
    p_S(S_obs == 0)  = 1; p_L(L_obs == 0)  = 1; p_N(N_obs == 0)  = 1;

    s            = struct();
    s.strain     = strain;
    s.TF         = TF_list;
    s.L_obs      = L_obs;
    s.N_obs      = N_obs;
    s.S_obs      = S_obs;
    s.p_S        = p_S;
    s.p_L        = p_L;
    s.p_N        = p_N;
    s.sig_S      = p_S < cfg.alpha;
    s.sig_L      = p_L < cfg.alpha;
    s.sig_N      = p_N < cfg.alpha;
    s.S_null_mean= sum_S / cfg.n_perm;
    s.TF_present = TF_present;
    s.sig_TFs    = TF_list(any(s.sig_S, 2));
    s.n_perm     = cfg.n_perm;
    s.alpha      = cfg.alpha;
    perm{k}      = s;

    fprintf('  %-8s significant TFs (p < %.2f at >= 1 transition): %d / %d  [%.1f s]\n\n', ...
        strain, cfg.alpha, numel(s.sig_TFs), nTF, toc(t_start));
end
end
