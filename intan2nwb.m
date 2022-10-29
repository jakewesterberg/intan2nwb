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
% Jake's github page (westerberg-science). Kilosort (2 is used here) is
% required for spike sorting, extended-GLM-for-synapse-detection (a modified
% version found on Jake's github) is required to estimate connectivity
% between units. Also, this of course requires the matnwb toolbox.

% Notes
% 1. This version of the code requires having a google sheet with some
% information pertaining to the recordings.

% NEED TO WRITE BETTER AMMENDMENT/RECALCULATION COMMANDS

function intan2nwb(varargin)
%% Defaults
skip_completed                  = false;
this_ident                      = []; % used to specify specific session(s) with their ident

%% pathing...can change to varargin or change function defaults for own machine
pp = pipelinePaths();

%% Varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'skip'}
            skip_completed = varargin{varStrInd(iv)+1};
        case {'this_ident'}
            this_ident = varargin{varStrInd(iv)+1};
        case {'ID'}
            ID = varargin{varStrInd(iv)+1};
    end
end

%% Use UI if desired
if ~exist('ID', 'var')
    ID = load(uigetfile(pwd, 'SELECT RECORDING ID FILE'));
end

%% Read recording session information
url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s', ID);
recording_info = webread(url_name);

% Create default processing list
to_proc = 1:length(unique(recording_info.Identifier));

% Limit to a specific session in a specific subject (ie identifier)
if ~isempty(this_ident)
    for ii = 1 : numel(this_ident)
        to_proc = find(strcmp(recording_info.Identifier, this_ident{ii}), 1);
    end
end

n_procd = 0;
%% Loop through sessions
for ii = to_proc

    % Skip files already processed if desired
    if exist([pp.NWB_DATA recording_info.Identifier{ii} '.nwb'], 'file') & skip_completed
        continue;
    end

    % Initialize nwb file
    nwb                                 = NwbFile;
    nwb.identifier                      = recording_info.Identifier{ii};
    nwb.session_start_time              = datetime(recording_info.Session(ii));
    nwb.general_experimenter            = recording_info.Investigator{ii};
    nwb.general_institution             = recording_info.Institution{ii};
    nwb.general_lab                     = recording_info.Lab{ii};
    nwb.general_session_id              = recording_info.Identifier{ii};
    nwb.general_experiment_description  = recording_info.Experiment_Description{ii};

    num_recording_devices = sum(strcmp(recording_info.Identifier, recording_info.Identifier{ii}));
    
    [nwb, recdev, probe] = i2nPRO(pp, nwb, recording_info, ii, num_recording_devices);

    probe_ctr = 0;
    for rd = 1 : num_recording_devices

        % Record analog traces
        nwb = i2nAIO(pp, nwb, recdev{rd});
        % Digital events
        nwb = i2nDIO(pp, nwb, recdev{rd});

        % Loop through probes to setup nwb tables
        for jj = 1 : recording_info.Probe_Count

            % BIN DATA
            if ~exist([pp.BIN_DATA nwb.identifier filesep ...
                    nwb.identifier '_probe-' num2str(probe{probe_ctr+1}.num) '.bin'], 'file')
                i2nBIN(pp, nwb, recdev{rd}, probe{probe_ctr+1});
            end

            % LFP AND MUA CALC
            try 
                eval(['nwb.acquisition.probe_' num2str(probe{probe_ctr+1}.num) '_lfp']);
            catch
                nwb = i2nCDS(pp, nwb, recdev{rd}, probe{probe_ctr+1});
            end

            % SPIKE SORTING
            if ~exist([pp.SPK_DATA nwb.identifier filesep ...
                    nwb.identifier filesep 'probe-' num2str(probe{probe_ctr+1}.num) ...
                    filesep 'rez2.mat'], 'file')
                nwb = i2nSPK(pp, nwb, recdev{rd}, probe{probe_ctr+1});
            end

            probe_ctr = probe_ctr + 1;

        end
    end
    n_procd = n_procd + 1;
end
disp(['SUCCESSFULLY PROCESSED ' n_procd ' FILES.'])
end