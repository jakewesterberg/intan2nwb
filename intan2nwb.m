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

function intan2nwb(ID, varargin)
%% Defaults - need to ammend
keepers                         = {};
skip_completed                  = true;
this_ident                      = []; % used to specify specific session(s) with their ident

%% pathing...can change to varargin or change function defaults for own machine
pp = pipelinePaths();
                                        
% add toolboxes
addpath(genpath(pp.TBOXES));

%% Varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'skip'}
            skip_completed = varargin{varStrInd(iv)+1};
        case {'this_ident'}
            this_ident = varargin{varStrInd(iv)+1};
        case {'SLACK_ID'}
            SLACK_ID = varargin{varStrInd(iv)+1};
    end
end

%% Prepare slack
if exist('SLACK_ID', 'var')
    send_slack_alerts = true;
else
    send_slack_alerts = false;
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
    if (exist([pp.DATA_DEST '_6_NWB_DATA' filesep recording_info.Identifier{ii} '.nwb'], 'file') & skip_completed) | ...
            recording_info.Preprocessor_Ignore_Flag(ii) == 1
        continue;
    end

    tic
    ttt=toc;
    if send_slack_alerts
        SendSlackNotification( ...
                SLACK_ID, ...
                [recording_info.Identifier{ii} ': [' s2HMS(ttt) '] Session processing started.'], ...
                'preprocess', ...
                'iJakebot', ...
                '', ...
                ':robot_face:');
    end

    if strcmp(recording_info.Raw_Data_Format{ii}, 'AI-NWB')

        if send_slack_alerts
            SendSlackNotification( ...
                SLACK_ID, ...
                [recording_info.Identifier{ii} ': [' s2HMS(ttt) '] Allen Institute NWB-raw-data identified.'], ...
                'preprocess', ...
                'iJakebot', ...
                '', ...
                ':robot_face:');
        end

        raw_data_dir = [pp.RAW_DATA 'dandi' filesep '000253' ...
            filesep 'sub_' recording_info.Subject{ii} ...
            filesep 'sub_' recording_info.Subject{ii} ...
            'sess_' num2str(recording_info.Session(ii)) ...
            filesep 'sub_' recording_info.Subject{ii} ...
            '+sess_' num2str(recording_info.Session(ii)) ...
            '_ecephys.nwb'];

        fpath = fileparts(raw_data_dir);

        if ~exist(raw_data_dir, 'file')
            warning('raw allen data not detected.')
            continue
        end

        copyfile(raw_data_dir, [pp.NWB_DATA recording_info.Identifier{ii} '.nwb'])
        nwb = nwbRead([pp.NWB_DATA recording_info.Identifier{ii} '.nwb']);

        i2nAIC(pp, nwb, recording_info, ii, fpath);

        if send_slack_alerts
            SendSlackNotification( ...
                SLACK_ID, ...
                [recording_info.Identifier{ii} ': [' s2HMS(ttt) '] AI conversion complete.'], ...
                'preprocess', ...
                'iJakebot', ...
                '', ...
                ':robot_face:');
        end

        n_procd = n_procd + 1;
        continue

    end

    % Initialize nwb file
    nwb                                 = NwbFile;
    nwb.identifier                      = recording_info.Identifier{ii};
    nwb.session_start_time              = datetime(datestr(datenum(num2str(recording_info.Session(ii)), 'yymmdd')));
    nwb.general_experimenter            = recording_info.Investigator{ii};
    nwb.general_institution             = recording_info.Institution{ii};
    nwb.general_lab                     = recording_info.Lab{ii};
    nwb.general_session_id              = recording_info.Identifier{ii};
    nwb.general_experiment_description  = recording_info.Experiment_Description{ii};

    num_recording_devices = sum(strcmp(recording_info.Identifier, recording_info.Identifier{ii}));

    for rd = 1 : num_recording_devices

        % RAW DATA
        fd1 = findDir(pp.RAW_DATA, nwb.identifier);
        fd2 = findDir(pp.RAW_DATA, ['dev-' num2str(rd-1)]);
        raw_data_present = sum(ismember(fd2, fd1));
        clear fd*

        if ~raw_data_present

            if (exist([pp.SCRATCH '\i2n_grab_data.bat'],'file'))
                delete([pp.SCRATCH '\i2n_grab_data.bat']);
            end

            fd1 = findDir(pp.DATA_SOURCE, nwb.identifier);
            fd2 = findDir(pp.DATA_SOURCE, ['dev-' num2str(rd-1)]);
            raw_data_temp = fd2(ismember(fd2, fd1));
            raw_data_temp = raw_data_temp{1};

            [~, dir_name_temp] = fileparts(raw_data_temp);

            % Grab data if missing
            workers = feature('numcores');
            fid = fopen([pp.SCRATCH '\i2n_grab_data.bat'], 'w');

            fprintf(fid, '%s\n', ...
                ['robocopy ' ...
                raw_data_temp ...
                ' ' ...
                [pp.RAW_DATA dir_name_temp] ...
                ' /e /j /mt:' ...
                num2str(workers)]);

            fclose('all');
            system([pp.SCRATCH '\i2n_grab_data.bat']);
            delete([pp.SCRATCH '\i2n_grab_data.bat']);
            ttt = toc;

            if send_slack_alerts
                SendSlackNotification( ...
                    SLACK_ID, ...
                    [nwb.identifier ': [' s2HMS(ttt) '] dev-' num2str(rd-1) ', Downloaded raw data from server.'], ...
                    'preprocess', ...
                    'iJakebot', ...
                    '', ...
                    ':robot_face:');
            end

        end
    end

    [nwb, recdev, probe] = i2nPRO(pp, nwb, recording_info, ii, num_recording_devices);
    ttt=toc;
    if send_slack_alerts
        SendSlackNotification( ...
                SLACK_ID, ...
                [nwb.identifier ': [' s2HMS(ttt) '] NWB initialization complete.'], ...
                'preprocess', ...
                'iJakebot', ...
                '', ...
                ':robot_face:');
    end
    
    probe_ctr = 0;
    for rd = 1 : num_recording_devices

        % Record analog traces
        nwb = i2nAIO(nwb, recdev{rd});
        ttt = toc;
        if send_slack_alerts
            SendSlackNotification( ...
                SLACK_ID, ...
                [nwb.identifier ': [' s2HMS(ttt) '] dev-' num2str(rd-1) ', Extracted analog I/O.'], ...
                'preprocess', ...
                'iJakebot', ...
                '', ...
                ':robot_face:');
        end

        % Digital events
        nwb = i2nDIO(nwb, recdev{rd});
        ttt = toc;
        if send_slack_alerts
            SendSlackNotification( ...
                SLACK_ID, ...
                [nwb.identifier ': [' s2HMS(ttt) '] dev-' num2str(rd-1) ', Extracted Digital I/O.'], ...
                'preprocess', ...
                'iJakebot', ...
                '', ...
                ':robot_face:');
        end

        % Loop through probes to setup nwb tables
        for jj = 1 : recording_info.Probe_Count(ii)

            % BIN DATA
            if ~exist([pp.BIN_DATA nwb.identifier filesep ...
                    nwb.identifier '_probe-' num2str(probe{probe_ctr+1}.num) '.bin'], 'file')
                i2nBIN(pp, nwb, recdev{rd}, probe{probe_ctr+1});
                ttt = toc;

                if send_slack_alerts
                    SendSlackNotification( ...
                        SLACK_ID, ...
                        [nwb.identifier ': [' s2HMS(ttt) '] dev-' num2str(rd-1) ', Binarized raw data.'], ...
                        'preprocess', ...
                        'iJakebot', ...
                        '', ...
                        ':robot_face:');
                end
            end

            % LFP AND MUA CALC
            try
                eval(['nwb.acquisition.probe_' num2str(probe{probe_ctr+1}.num) '_lfp']);
            catch
                nwb = i2nCDS(nwb, recdev{rd}, probe{probe_ctr+1});
                ttt = toc;

                if send_slack_alerts
                    SendSlackNotification( ...
                        SLACK_ID, ...
                        [nwb.identifier ': [' s2HMS(ttt) '] dev-' num2str(rd-1) ', Filtered MUAe and LFP.'], ...
                        'preprocess', ...
                        'iJakebot', ...
                        '', ...
                        ':robot_face:');
                end

            end

            % SPIKE SORTING
            try
                nwb = i2nSPK(pp, nwb, recdev{rd}, probe{probe_ctr+1});

                ttt = toc;

                if send_slack_alerts
                    SendSlackNotification( ...
                        SLACK_ID, ...
                        [nwb.identifier ': [' s2HMS(ttt) '] dev-' num2str(rd-1) ', Completed spike sorting/curation (~' ...
                        num2str(sum(nwb.units.vectordata.get('quality').data(:))) ' good units).'], ...
                        'preprocess', ...
                        'iJakebot', ...
                        '', ...
                        ':robot_face:');
                end
            catch
                warning('KS DIDNT PAN OUT!!!!!!')
            end

            probe_ctr = probe_ctr + 1;

        end
    end

    ttt = toc;
    if send_slack_alerts
        SendSlackNotification( ...
            SLACK_ID, ...
            [nwb.identifier ': [' s2HMS(ttt) '] Session processing complete.'], ...
            'preprocess', ...
            'iJakebot', ...
            '', ...
            ':robot_face:');
    end

    n_procd = n_procd + 1;
    nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

    % Cleanup
    i2nCleanup(pp, keepers);


end
disp(['SUCCESSFULLY PROCESSED ' num2str(n_procd) ' FILES.'])
end