function TF_target_info = build_tf_target_info(data_info, tflink, cfg)
%BUILD_TF_TARGET_INFO  Strain-specific TF universe and TFLink sub-network.
%
%   TF_target_info = BUILD_TF_TARGET_INFO(data_info, tflink, cfg)
%
%   TF_target_info{k,1} : [nTF x 1]   string. Union over all timepoints of the
%                         transcription factors detected in strain k. This is
%                         the set of TFs that can ever act as a source node.
%   TF_target_info{k,2} : [nEdge x 2] string. The TFLink edges whose source is
%                         one of those TFs, i.e. the time-invariant regulatory
%                         backbone available to strain k.
%
%   Identifiers are normalized (see NORMALIZE_IDS) so that downstream set
%   operations are case-insensitive.

TF_target_info = cell(cfg.n_strain, 2);

src = normalize_ids(tflink(:, 1));
dst = normalize_ids(tflink(:, 2));

for k = 1:cfg.n_strain
    tf_union = strings(0, 1);
    for t = 1:cfg.n_tp
        tf_t     = normalize_ids(data_info{t, 2}{k, cfg.col_tf}(:, 1));
        tf_union = union(tf_union, tf_t);
    end
    tf_union(tf_union == "") = [];
    TF_target_info{k, 1} = tf_union;

    keep = ismember(src, tf_union);
    TF_target_info{k, 2} = [src(keep), dst(keep)];

    if cfg.verbose
        fprintf('%-8s : %5d TFs detected, %7d TFLink edges retained\n', ...
            cfg.strain_labels{k}, numel(tf_union), sum(keep));
    end
end
if cfg.verbose; fprintf('\n'); end
end
