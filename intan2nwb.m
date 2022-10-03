%% Header
% Jake Westerberg, PhD (westerberg-science)
% Vanderbilt University
% jakewesterberg@gmail.com
% Code contributions from Patrick Meng (VU)

% Description
% Written as a pipeline for data collected in the Bastos Lab (or similar
% setup using the intan system) to be transformed from raw to the nwb
% format.

% Requirements
% Certain aspects of the data processing require toolboxes found in other
% github repos. Original or forked versions of all required can be found on
% Jake's github page (westerberg-science). Kilosort (3 is used here) is
% required for spike sorting, potential-toolbox is required for most low
% frequency computations, extended-GLM-for-synapse-detection (a modified
% version found on Jake's github) is required to estimate connectivity
% between units. Also, this of course requires the matnwb toolbox.

% Notes
% 1. This version of the code requires having a google sheet with some
% information pertaining to the recordings. 


function intan2nwb(ID, varargin)
%% Defaults
gpu_compute                     = true; % use gpu processing where possible
parallel_compute                = true; % use parallel computing where possible

skip_completed                  = true;

in_file_path                    = '\\teba.psy.vanderbilt.edu\bastoslab\_BL_DATA_PIPELINE\_RAW_DATA\';
out_file_path                   = '\\teba.psy.vanderbilt.edu\bastoslab\_BL_DATA_PIPELINE\_NWB_DATA\';

this_subject                    = []; % used to specify processing for only certain subjects
this_ident                      = []; % used to specify specific session(s) with their ident

params.save_bin                 = true;

params.spk                      = true; % requires Kilosort (v3) on path currently only option
params.cnx                      = true; % requires extended-GLM-for-synapse-detection on path
params.muae                     = true; % requires potential-toolbox on path
params.lfp                      = true;
params.csd                      = true; % requires potential-toolbox on path
params.nflfp                    = true; % requires potential-toolbox on path
params.eye                      = true;

reprocess.spk                   = true;
reprocess.cnx                   = true;

params.trigger_PRO              = false; % not useable now.

%% Varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-i', 'in_file_path'}
            out_file_path = varargin{varStrInd(iv)+1};
        case {'-o', 'out_file_path'}
            out_file_path = varargin{varStrInd(iv)+1};
        case {'skip'}
            skip_completed = varargin{varStrInd(iv)+1};
        case {'this_subject'}
            this_subject = varargin{varStrInd(iv)+1};      
        case {'this_ident'}
            this_ident = varargin{varStrInd(iv)+1};
        case {'-p', 'params'}
            params = varargin{varStrInd(iv)+1};
        case {'-pc', 'parallel_compute'}
            parallel_compute = varargin{varStrInd(iv)+1};
        case {'-gpu', 'gpu_compute'}
            gpu_compute = varargin{varStrInd(iv)+1};          
    end
end

%% Read recording session information
url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s', ID);
recording_info = webread(url_name);

% Create default processing list
n_idents = length(recording_info.Identifier);
to_proc = 1:n_idents;

% Limit to sessions within subject (if applicable)
if ~isempty(this_subject)
    to_proc = find(strcmp(recording_info.Subject, this_subject));
end

% Limit to a specific session in a specific subject
if ~isempty(this_ident)
    to_proc = nan(1, numel(this_ident));
    for ii = 1 : numel(this_ident)
        to_proc(ii) = find(strcmp(recording_info.Identifier, this_ident{ii}));
    end
end

%% Loop through sessions
for ii = to_proc
    
    % Skip files already processed if desired
    if exist([out_file_path 'sub-' recording_info.Subject '_ses-' recording_info.Session '.nwb'], 'file') & ...
            skip_completed
        continue;
    end

    % Initialize nwb file
    nwb = NwbFile;

    % Read settings
    settings_file = find_in_dir(in_file_path, 'settings.xml');
    try; settings_info = parseXML(settings_file{1});
    catch; error('intan recording does not have readable settings file in directory');
    end

    % Record subject information
    nwb.identifier = recording_info.Identifier(ii);

    % Record spiking data
    if params.spk
        if exist('', 'file') & ~reprocess.spk
        else
        end
    end

    % Estimate interconnectivity of single units
    if params.cnx

    end

    % Record multiunit envelope
    if params.muae

    end

    % Record LFP
    if params.lfp

    end

    % Compute CSD
    if params.csd

    end

    % Compute nfLFP
    if params.nflfp

    end

    % Record eye trace
    if params.eye

    end

    % Save to NWB
    nwbExport(nwb, [out_file_path 'sub-' recording_info.Subject '_ses-' recording_info.Session '.nwb']);
    disp(['SUCCESSFULLY SAVED: ' out_file_path 'sub-' recording_info.Subject '_ses-' recording_info.Session '.nwb'])

    % Increment counter
    n_procd = n_procd + 1;

end

disp(['SUCCESSFULLY PROCESSED ' n_procd ' FILES.'])

end