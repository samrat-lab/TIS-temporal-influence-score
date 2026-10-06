function cfg = tis_config(varargin)
%TIS_CONFIG  Central configuration for the TIS pipeline.
%
%   cfg = TIS_CONFIG()
%   cfg = TIS_CONFIG('n_perm', 200, 'run_permutation', false)
%
%   Any field of the default configuration can be overridden with a
%   name-value pair. The defaults below reproduce the analysis reported in
%   the manuscript.

here       = fileparts(mfilename('fullpath'));
cfg        = struct();
cfg.root   = fileparts(here);

% ---- paths -------------------------------------------------------------
% data_info_file and tflink_file follow data_dir unless they are overridden
% explicitly, so either of these works:
%     run_all('data_dir', 'E:\...\New')
%     run_all('data_info_file', 'E:\...\New\data_info.mat', ...
%             'tflink_file',    'E:\...\Analysis_Matlab\Final_TF_Target.mat')
cfg.data_dir        = fullfile(cfg.root, 'data');
cfg.results_dir     = fullfile(cfg.root, 'results');
cfg.data_info_file  = '';    % resolved below
cfg.tflink_file     = '';    % resolved below

% ---- experimental design ----------------------------------------------
% Row order of data_info{t,2} must match strain_labels.
cfg.strain_labels = {'H37Ra', 'H37Rv', 'BND433', 'JAL2287'};
cfg.tp_hours      = [6 10 14 18 22 26 30 34 38 42 46 50];   % hours post-infection

% ---- column map of data_info{t,2} (one row per strain) -----------------
% 1 : master protein list
% 2 : newly synthesized proteins detected at this timepoint in this strain
% 3 : up-regulated proteins   (vs uninfected control)
% 4 : down-regulated proteins (vs uninfected control)
% 5 : DENSP = union of columns 3 and 4
% 6 : transcription factors detected at this timepoint
cfg.col_master       = 1;
cfg.col_newly_synth  = 2;
cfg.col_up           = 3;
cfg.col_down         = 4;
cfg.col_densp        = 5;
cfg.col_tf           = 6;

% The TIS definition uses the newly synthesized proteome (column 2) as the
% node set P_i. Set to cfg.col_densp to score on DENSP instead.
cfg.col_protein_pool = cfg.col_newly_synth;

% ---- permutation test ---------------------------------------------------
cfg.run_permutation = true;
cfg.n_perm          = 1000;    % 1000 in the manuscript; lower for quick runs
cfg.alpha           = 0.05;
cfg.rng_seed        = 42;

% ---- misc ---------------------------------------------------------------
cfg.reach_chunk = 500;   % rows per distances() block (memory control)
cfg.verbose     = true;

% ---- apply overrides ----------------------------------------------------
if mod(numel(varargin), 2) ~= 0
    error('tis_config:badArgs', 'Options must be given as name-value pairs.');
end
for a = 1:2:numel(varargin)
    name = varargin{a};
    if ~isfield(cfg, name)
        error('tis_config:unknownOption', 'Unknown option "%s".', name);
    end
    cfg.(name) = varargin{a+1};
end

% ---- resolve input paths ------------------------------------------------
if isempty(cfg.data_info_file)
    cfg.data_info_file = fullfile(cfg.data_dir, 'data_info.mat');
end
if isempty(cfg.tflink_file)
    cfg.tflink_file = fullfile(cfg.data_dir, 'Final_TF_Target.mat');
end

cfg.n_tp     = numel(cfg.tp_hours);
cfg.n_strain = numel(cfg.strain_labels);
end
