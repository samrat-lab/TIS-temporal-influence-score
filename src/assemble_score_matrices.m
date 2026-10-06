function scores = assemble_score_matrices(tis, TF_target_info, cfg)
%ASSEMBLE_SCORE_MATRICES  TF x timepoint matrices of the 1st and 2nd TIS terms.
%
%   scores = ASSEMBLE_SCORE_MATRICES(tis, TF_target_info, cfg)
%
%   Reshapes the per-timepoint score vectors returned by COMPUTE_TIS_SCORES
%   into one matrix per strain: rows are the strain TF universe, columns are
%   the timepoints.
%
%   The matrices are zero-initialised. A (TF, timepoint) pair scores 0 either
%   because the TF was not detected at that timepoint, or because none of its
%   TFLink targets were measured at the next timepoint - in both cases
%   L = N = 0 follows from the definition, so zero is the correct value and
%   not a missing one.
%
%   Returned structure, per strain k:
%     scores{k}.strain char
%     scores{k}.TF     [nTF x 1]    string - strain TF universe
%     scores{k}.L      [nTF x n_tp] double - 1st term (local influence)
%     scores{k}.N      [nTF x n_tp] double - 2nd term (propagated influence)
%     scores{k}.S      [nTF x n_tp] double - total TIS = L + N
%
%   The last column is always zero: no transition starts at the final
%   timepoint, so no score is defined there.

scores = cell(cfg.n_strain, 1);
nT     = cfg.n_tp;

for k = 1:cfg.n_strain
    TF  = TF_target_info{k, 1};
    nTF = numel(TF);

    L = zeros(nTF, nT);
    N = zeros(nTF, nT);

    for t = 1:nT
        tf_t = tis{k}.TF{t};
        if isempty(tf_t); continue; end

        [found, row] = ismember(tf_t, TF);
        assert(all(found), ['assemble_score_matrices: the TF list at ' ...
            'timepoint %d is not a subset of the strain TF universe.'], t);

        L(row, t) = tis{k}.L{t};
        N(row, t) = tis{k}.N{t};
    end

    s        = struct();
    s.strain = cfg.strain_labels{k};
    s.TF     = TF;
    s.L      = L;
    s.N      = N;
    s.S      = L + N;
    scores{k}= s;

    if cfg.verbose
        fprintf('%-8s : score matrix %dx%d, %d entries with S > 0\n', ...
            cfg.strain_labels{k}, nTF, nT, sum(s.S(:) > 0));
    end
end
if cfg.verbose; fprintf('\n'); end
end
