function intan2nwb(ID, varargin)
%% Header info
% Jake Westerberg (westerberg-science)
% Vanderbilt University
% jakewesterberg@gmail.com
% Code contributions from Patrick Meng

%% Defaults
skip_completed                  = true;

in_file_path                    = '\\teba.psy.vanderbilt.edu\bastoslab\_RAW_DATA\';
out_file_path                   = '\\teba.psy.vanderbilt.edu\bastoslab\_NWB_DATA\';

this_subject                    = [];
this_session                    = [];

params.binary_output            = true;
params.kilosort                 = true;
params.estimate_connectivity    = true;
params.muae                     = true;
params.lfp                      = true;
params.csd                      = true;
params.nflfp                    = true;
params.eye_tracking             = true;

%% Varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-i', 'in_file_path'}
            out_file_path = varargin{varStrInd(iv)+1};
        case {'-o', 'out_file_path'}
            out_file_path = varargin{varStrInd(iv)+1};
        case {'-a', 'subj'}
            subj = varargin{varStrInd(iv)+1};
        case {'-s', 'ses'}
            ses = varargin{varStrInd(iv)+1};
        case {'skip'}
            skip_completed = varargin{varStrInd(iv)+1};
        case {'this_subject'}
            this_subject = varargin{varStrInd(iv)+1};      
        case {'this_session'}
            this_session = varargin{varStrInd(iv)+1};
        case {'-p', 'params'}
            params = varargin{varStrInd(iv)+1};   
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
if ~isempty(this_session)
    to_proc = find(find(strcmp(recording_info.Subject, this_subject)) & ...
        find(strcmp(recording_info.Session, this_session)));
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
    if params.kilosort

    end

    % Estimate interconnectivity of single units
    if params.estimate_connectivity

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
    if params.eye_tracking

    end

    % Save to NWB
    nwbExport(nwb, [out_file_path 'sub-' recording_info.Subject '_ses-' recording_info.Session '.nwb']);

    % Increment counter
    n_procd = n_procd + 1;

end

disp(['SUCCESSFULLY PROCESSED ' n_procd ' FILES.'])

end