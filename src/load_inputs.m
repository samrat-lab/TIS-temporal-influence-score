function [data_info, tflink] = load_inputs(cfg)
%LOAD_INPUTS  Load data_info.mat and the TFLink edge list, with validation.
%
%   [data_info, tflink] = LOAD_INPUTS(cfg)
%
%   data_info : {n_tp x 2} cell. data_info{t,2} is a {n_strain x >=6} cell
%               of string columns (see data/README.md).
%   tflink    : [n_edge x 2] string. Column 1 = source TF, column 2 = target.

hint = ['\nEither copy the file into the data/ folder, or point the ' ...
        'pipeline at its\ncurrent location, e.g.\n' ...
        '    run_all(''data_dir'', ''C:\\path\\to\\folder'')\n' ...
        '    run_all(''data_info_file'', ''C:\\path\\to\\data_info.mat'', ...\n' ...
        '            ''tflink_file'', ''C:\\path\\to\\Final_TF_Target.mat'')\n' ...
        'See data/README.md for the expected contents.\n'];

if ~isfile(cfg.data_info_file)
    error('load_inputs:missingFile', ...
        ['data_info.mat not found at:\n  %s\n' hint], cfg.data_info_file);
end
if ~isfile(cfg.tflink_file)
    error('load_inputs:missingFile', ...
        ['Final_TF_Target.mat not found at:\n  %s\n' hint], cfg.tflink_file);
end

S = load(cfg.data_info_file);
if ~isfield(S, 'data_info')
    error('load_inputs:badVariable', ...
        '%s must contain a variable named "data_info".', cfg.data_info_file);
end
data_info = S.data_info;

T  = load(cfg.tflink_file);
if isfield(T, 'Final_TF_Target')
    tflink = T.Final_TF_Target;
else
    fn = fieldnames(T);
    if numel(fn) ~= 1
        error('load_inputs:badVariable', ...
            '%s must contain a variable named "Final_TF_Target".', cfg.tflink_file);
    end
    tflink = T.(fn{1});
end
tflink = string(tflink);

% ---- validation ---------------------------------------------------------
if size(data_info, 1) < cfg.n_tp
    error('load_inputs:badShape', ...
        'data_info has %d rows but %d timepoints are configured.', ...
        size(data_info, 1), cfg.n_tp);
end
if size(tflink, 2) < 2
    error('load_inputs:badShape', 'Final_TF_Target must have >= 2 columns.');
end
tflink = tflink(:, 1:2);

for t = 1:cfg.n_tp
    blk = data_info{t, 2};
    if size(blk, 1) < cfg.n_strain || size(blk, 2) < cfg.col_tf
        error('load_inputs:badShape', ...
            ['data_info{%d,2} is %dx%d but at least %dx%d is required ' ...
             '(strains x columns).'], t, size(blk,1), size(blk,2), ...
             cfg.n_strain, cfg.col_tf);
    end
end

if cfg.verbose
    fprintf('Loaded %d timepoints x %d strains; TFLink edges: %d\n', ...
        cfg.n_tp, cfg.n_strain, size(tflink, 1));
end
end
