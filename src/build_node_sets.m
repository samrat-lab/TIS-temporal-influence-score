function nodes = build_node_sets(data_info, cfg)
%BUILD_NODE_SETS  Per-timepoint protein and TF node sets for each strain.
%
%   nodes = BUILD_NODE_SETS(data_info, cfg)
%
%   nodes{k}.prot{t}  : V_t  - all proteins measured at timepoint t
%                              (newly synthesized proteome, data_info column 2)
%   nodes{k}.tf{t}    : T_t  - transcription factors detected at timepoint t
%                              (data_info column 6)
%   nodes{k}.nonTF{t} : P_t  = V_t \ T_t - newly synthesized non-TF proteins
%
%   All sets are unique and normalized. P_t is the denominator of the local
%   influence term L (Eq. 10).

nodes = cell(cfg.n_strain, 1);

for k = 1:cfg.n_strain
    s       = struct();
    s.prot  = cell(cfg.n_tp, 1);
    s.tf    = cell(cfg.n_tp, 1);
    s.nonTF = cell(cfg.n_tp, 1);

    for t = 1:cfg.n_tp
        prot = normalize_ids(data_info{t, 2}{k, cfg.col_protein_pool}(:, 1));
        tf   = normalize_ids(data_info{t, 2}{k, cfg.col_tf}(:, 1));
        prot(prot == "") = [];
        tf(tf   == "")   = [];

        s.prot{t}  = unique(prot);
        s.tf{t}    = unique(tf);
        s.nonTF{t} = setdiff(s.prot{t}, s.tf{t});
    end
    nodes{k} = s;

    if cfg.verbose
        fprintf('%-8s : proteins/tp = [%s]\n', cfg.strain_labels{k}, ...
            strjoin(string(cellfun(@numel, s.prot)'), ' '));
        fprintf('%-8s : TFs/tp      = [%s]\n', '', ...
            strjoin(string(cellfun(@numel, s.tf)'), ' '));
    end
end
if cfg.verbose; fprintf('\n'); end
end
