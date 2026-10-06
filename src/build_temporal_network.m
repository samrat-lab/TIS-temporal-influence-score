function net = build_temporal_network(nodes, TF_target_info, cfg)
%BUILD_TEMPORAL_NETWORK  Directed, layered TF -> target network (Eqs. 1-4).
%
%   net = BUILD_TEMPORAL_NETWORK(nodes, TF_target_info, cfg)
%
%   Edges run strictly forward in time, from an ACTIVE transcription factor at
%   layer i to one of its TFLink targets measured at layer i+1:
%
%       E_i = { (t,p) : t in A_i, p in V_{i+1}, (t -> p) in TFLink }
%
%   The active TF set is propagated through the layers:
%
%       A_1     = TFs detected at t1 that are TFLink sources
%       A_{i}   = { TF targets reached at layer i that are TFLink sources }
%                 union
%                 { TFs detected at layer i that are TFLink sources }
%
%   so that a TF enters the network either because it was reached from the
%   previous layer or because it was newly detected at that layer.
%
%   Returned structure, per strain k:
%     net{k}.edges      [E x 2] string - tagged node ids "NAME_<layer>"
%     net{k}.active{i}  [a_i x 1] string - active TF set at layer i (untagged)
%     net{k}.layer_edges{i} - edges of transition i -> i+1 (tagged)

net = cell(cfg.n_strain, 1);
nT  = cfg.n_tp;

for k = 1:cfg.n_strain
    src        = TF_target_info{k, 2}(:, 1);
    dst        = TF_target_info{k, 2}(:, 2);
    tf_sources = unique(src);

    active      = cell(nT, 1);
    layer_edges = cell(nT - 1, 1);

    active{1} = intersect(nodes{k}.tf{1}, tf_sources);

    for i = 2:nT
        mask = ismember(src, active{i-1});
        es   = src(mask);
        ed   = dst(mask);

        keep = ismember(ed, nodes{k}.prot{i});       % target measured at layer i
        es   = es(keep);
        ed   = ed(keep);

        e = unique([es, ed], 'rows');
        if isempty(e)
            layer_edges{i-1} = strings(0, 2);
        else
            layer_edges{i-1} = [tag_nodes(e(:,1), i-1), tag_nodes(e(:,2), i)];
        end

        reached_tf  = intersect(unique(ed), tf_sources);
        detected_tf = intersect(nodes{k}.tf{i}, tf_sources);
        active{i}   = union(reached_tf, detected_tf);
    end

    s             = struct();
    s.active      = active;
    s.layer_edges = layer_edges;
    s.edges       = vertcat(layer_edges{:});
    if isempty(s.edges); s.edges = strings(0, 2); end
    net{k} = s;

    if cfg.verbose
        fprintf('%-8s : %d temporal edges, %d nodes\n', cfg.strain_labels{k}, ...
            size(s.edges, 1), numel(unique(s.edges(:))));
    end
end
if cfg.verbose; fprintf('\n'); end
end
