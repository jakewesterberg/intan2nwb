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


function intan2nwb(varargin)
%% Defaults
gpu_compute                     = true; %true; % use gpu processing where possible
workers                         = 1; % use parallel computing where possible

skip_completed                  = true;

in_file_path                    = '\\teba.psy.vanderbilt.edu\bastoslab\_BL_DATA_PIPELINE\_0_RAW_DATA\';
out_file_path                   = '\\teba.psy.vanderbilt.edu\bastoslab\_BL_DATA_PIPELINE\_2_NWB_DATA\';

this_subject                    = []; % used to specify processing for only certain subjects
this_ident                      = []; % used to specify specific session(s) with their ident

bin_file_path                    = '\\teba.psy.vanderbilt.edu\bastoslab\_BL_DATA_PIPELINE\_1_BIN_DATA\';
params.save_bin                 = true;

params.downsample_fs            = 1000;

params.spk                      = false; % requires Kilosort (v3) on path currently only option
params.cnx                      = false; % requires extended-GLM-for-synapse-detection on path
params.csd                      = false; % requires potential-toolbox on path
params.nflfp                    = false; % requires potential-toolbox on path
params.eye                      = false;

reprocess.spk                   = true;
reprocess.cnx                   = true;

params.trigger_PRO              = false; % not useable now.

%% Varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-i', 'in_file_path'}
            in_file_path = varargin{varStrInd(iv)+1};
        case {'-o', 'out_file_path'}
            out_file_path = varargin{varStrInd(iv)+1};
        case {'-b', 'bin_file_path'}
            bin_file_path = varargin{varStrInd(iv)+1};
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
        case {'ID'}
            ID = varargin{varStrInd(iv)+1};
    end
end

%% Use UI if desired
if ~exist('ID', 'var')
    ID = load(uigetfile(pwd, 'SELECT RECORDING ID FILE'));
    in_file_path = uigetdir(pwd, 'SELECT DATA INPUT DIRECTORY');
    out_file_path = uigetdir(pwd, 'SELECT DATA OUTPUT DIRECTORY');
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

    % Find the correct subpath
    in_file_path_1 = findDir(in_file_path, datestr(recording_info.Session(ii), 'yymmdd'));
    in_file_path_2 = findDir(in_file_path, recording_info.Subject{ii});

    match_dir = ismember(in_file_path_2, in_file_path_1);
    if isempty(match_dir)
        in_file_path_2 = findDir(in_file_path, recording_info.Subject_Nickname{ii});
        match_dir = ismember(in_file_path_2, in_file_path_1);
        if isempty(match_dir)
            warning(['COULD NOT FIND DIR FOR ' recording_info.Subject{ii} '-' ...
                datestr(recording_info.Session(ii), 'yymmdd') ' MOVING ON.'])
            continue
        end
    end

    in_file_path_itt = [in_file_path_2{match_dir} filesep];
    clear in_file_path_1 in_file_path_2 match_dir

    % Create file identifier
    file_ident = ['sub-' recording_info.Subject{ii} '_ses-' datestr(recording_info.Session(ii), 'yymmdd')];

    % Skip files already processed if desired
    if exist([out_file_path file_ident '.nwb'], 'file') & ...
            skip_completed
        continue;
    end

    % Read settings
    intan_header = readIntanHeader(in_file_path_itt);

    % Initialize nwb file
    nwb                                 = NwbFile;
    nwb.identifier                      = recording_info.Identifier{ii};
    nwb.session_start_time              = datetime(recording_info.Session(ii));
    nwb.general_experimenter            = recording_info.Investigator{ii};
    nwb.general_institution             = recording_info.Institution{ii};
    nwb.general_lab                     = recording_info.Lab{ii};
    nwb.general_session_id              = recording_info.Identifier{ii};
    nwb.general_experiment_description  = recording_info.Experiment_Description{ii};

    % Determine which probes are present
    probes = strtrim(split(recording_info.Probe_Ident{ii}, ','));

    % Initialize probe table
    variables = {'x', 'y', 'z', 'imp', 'location', 'filtering', 'group', 'label'};
    e_table = cell2table(cell(0, length(variables)), 'VariableNames', variables);

    % Loop through probes to setup nwb tables
    for jj = 1 : recording_info.Probe_Count

        % Determine number of channels that should be present for probe
        if strcmp(class(recording_info.Probe_Channels), 'double')
            n_channels = recording_info.Probe_Channels(jj);
        elseif strcmp(class(recording_info.Probe_Channels), 'cell')
            temp_array_1 = strtrim(split(recording_info.Probe_Channels{ii}, ','));
            n_channels = str2double(temp_array_1{jj});
            clear temp_array_1
        end

        % Load the correct channel map file
        load([probes{jj} '.mat'], 'channel_map', 'x', 'y')

        % Setup some basic NWB information
        if strcmp(class(recording_info.Probe_Shanks), 'double')
            n_shanks = recording_info.Probe_Shanks(jj);
        elseif strcmp(class(recording_info.Probe_Shanks), 'cell')
            temp_array_1 = strtrim(split(recording_info.Probe_Shanks{ii}, ','));
            n_shanks = str2double(temp_array_1{jj});
            clear temp_array_1
        end

        if strcmp(class(recording_info.Probe_Channels_Per_Shank), 'double')
            n_channels_per_shank = recording_info.Probe_Channels_Per_Shank(jj);
        elseif strcmp(class(recording_info.Probe_Channels_Per_Shank), 'cell')
            temp_array_1 = strtrim(split(recording_info.Probe_Channels_Per_Shank{ii}, ','));
            n_channels_per_shank = str2double(temp_array_1{jj});
            clear temp_array_1
        end

        % Create device
        device = types.core.Device(...
            'description', paren(strtrim(split(recording_info.Probe_Ident{ii}, ',')), jj), ...
            'manufacturer', paren(strtrim(split(recording_info.Probe_Manufacturer{ii}, ',')), jj) ...
            );

        % Input device information
        nwb.general_devices.set('array', device);
        for ishank = 1:n_shanks
            electrode_group = types.core.ElectrodeGroup( ...
                'description', ['electrode group for shank' num2str(ishank)], ...
                'location', paren(strtrim(split(recording_info.Area{ii}, ',')), jj), ...
                'device', types.untyped.SoftLink(device) ...
                );
            nwb.general_extracellular_ephys.set(['shank' num2str(ishank)], electrode_group);
            group_object_view = types.untyped.ObjectView(electrode_group);

            % Grab X position
            if strcmp(class(recording_info.X), 'double')
                temp_X = recording_info.X(jj);
            elseif strcmp(class(recording_info.X), 'cell')
                temp_array_1 = strtrim(split(recording_info.X{ii}, ','));
                temp_X = str2double(temp_array_1{jj});
                clear temp_array_1
            end

            % Grab Y position
            if strcmp(class(recording_info.Y), 'double')
                temp_Y = recording_info.Y(jj);
            elseif strcmp(class(recording_info.X), 'cell')
                temp_array_1 = strtrim(split(recording_info.Y{ii}, ','));
                temp_Y = str2double(temp_array_1{jj});
                clear temp_array_1
            end

            % Get Z position initialized
            if strcmp(class(recording_info.Z), 'double')
                init_Z = recording_info.Z(jj);
            elseif strcmp(class(recording_info.Z), 'cell')
                temp_array_1 = strtrim(split(recording_info.Z{ii}, ','));
                init_Z = str2double(temp_array_1{jj});
                clear temp_array_1
            end

            temp_imp = NaN; % Can we add in impedance data, day-to-day from intan file?

            temp_loc = paren(strtrim(split(recording_info.Area{ii}, ',')), jj);

            temp_filt = NaN; % Can probably grab this from the settings file eventually.

            for ielec = 1:n_channels_per_shank
                electrode_label = ['shank' num2str(ishank) 'elec' num2str(ielec)];

                temp_Z = init_Z - (max(abs(y)) - y(ielec));

                e_table = [e_table; {temp_X, temp_Y, temp_Z, temp_imp, temp_loc, temp_filt, group_object_view, electrode_label}];
            end
        end
    end

    % Record electrode table
    electrode_table = util.table2nwb(e_table, 'all electrodes');
    nwb.general_extracellular_ephys_electrodes = electrode_table;

    % Initialize electrode table region
    electrode_table_region = types.hdmf_common.DynamicTableRegion( ...
        'table', types.untyped.ObjectView(electrode_table), ...
        'description', 'all electrodes', ...
        'data', (0:height(e_table)-1)');

    % Determine number of samples in datafiles
    n_samples = length(intan_header.time_stamp);

    % Determine the downsampling
    params.downsample_factor = intan_header.sampling_rate/params.downsample_fs;
    downsample_size = length(downsample(intan_header.time_stamp, params.downsample_factor));

    % Initialize DC offset filter
    [DC_offset_bwb, DC_offset_bwa] = butter(1, 0.1/(intan_header.sampling_rate/2), 'high');

    % Initialize filter information
    [muae_bwb, muae_bwa] = butter(2, [500 5000]/(intan_header.sampling_rate/2), 'bandpass');
    [muae_power_bwb, muae_power_bwa] = butter(4, 250/(intan_header.sampling_rate/2), 'low');
    muae_all = [];

    [lfp_bwb, lfp_bwa] = butter(2, [1 250]/(intan_header.sampling_rate/2), 'bandpass');
    lfp_all = [];
    if params.csd || params.nflfp
        csd_all = [];
        if params.nflfp
            nflfp_all = [];
        end
    end

    % Loop through the individual probe data
    for jj = 1 : recording_info.Probe_Count

        % Determine number of channels that should be present for probe
        if strcmp(class(recording_info.Probe_Channels), 'double')
            n_channels = recording_info.Probe_Channels(jj);
        elseif strcmp(class(recording_info.Probe_Channels), 'cell')
            temp_array_1 = strtrim(split(recording_info.Probe_Channels{ii}, ','));
            n_channels = str2double(temp_array_1{jj});
            clear temp_array_1
        end

        % Load the correct channel map file
        load([probes{jj} '.mat'], 'channel_map', 'x', 'y')

        % Initialize data matrices. Need to fix for multiprobe
        lfp = zeros(n_channels, downsample_size);
        muae = zeros(n_channels, downsample_size);
        if params.csd;      csd = zeros(n_channels, downsample_size);    end
        if params.nflfp;    nflfp = zeros(n_channels, downsample_size);  end

        %pool1 = parpool(workers);
        for kk = 1:n_channels

            % Open file and init data
            current_fid             = fopen(in_file_path_itt + "\amp-" + intan_header.amplifier_channels(kk).native_channel_name + ".dat");

            % Setup array on GPU or in mem depending on run parameters
            if gpu_compute
                current_data            = gpuArray(double(fread(current_fid, n_samples, 'int16')) * 0.195);
            else
                current_data            = double(fread(current_fid, n_samples, 'int16')) * 0.195;
            end

            % Do DC offset filter
            current_data    = filtfilt(DC_offset_bwb, DC_offset_bwa, current_data);

            % Do data type specific filtering
            if gpu_compute % Yes, I know this could be written better
                muae(kk,:)  = gather(downsample(filtfilt(muae_power_bwb, muae_power_bwa, ...
                    abs(filtfilt(muae_bwb, muae_bwa, current_data))), params.downsample_factor));
                lfp(kk,:)   = gather(downsample(filtfilt(lfp_bwb, lfp_bwa, current_data), params.downsample_factor));


            else
                muae(kk,:)  = downsample(filtfilt(muae_power_bwb, muae_power_bwa, ...
                    abs(filtfilt(muae_bwb, muae_bwa, current_data))), params.downsample_factor);
                lfp(kk,:)   = downsample(filtfilt(lfp_bwb, lfp_bwa, current_data), params.downsample_factor);
            end

            % Close file
            fclose(current_fid);
            disp([num2str(kk) '/' num2str(n_channels) ' COMPLETED.'])
        end
        %delete(pool1)

        %Rearrange the channels to the order on the probe (starts at 0, +1 so it
        %matches matlab indexing)
        muae = muae(channel_map+1,:);
        muae_all = [muae_all; muae];

        lfp = lfp(channel_map+1,:);
        lfp_all = [lfp_all; lfp];


        if params.csd || params.nflfp

            % CSD computation currently uses Nicholson and Freeman
            % method. Will create a flag to use iCSD when that is
            % properly implemented.
            csd = CSD(lfp, y);

            if params.csd
                csd_all = [csd_all; csd];
            end

            if params.nflfp
                nflfp = NFLFP(csd, y);
                nflfp_all = [nflfp_all; nflfp];
            end
        end
    end

    % Record spiking data
    if params.spk
        if exist('', 'file') & ~reprocess.spk
        else
            for kk = unique(cat(1,intan_header.amplifier_channels.port_prefix))
                intan2bin(in_file_path_itt, bin_file_path, [file_ident '.bin'], kk)
            end
        end
    end

    % Initialize electrical series
    lfp_electrical_series = types.core.ElectricalSeries( ...
        'starting_time', 0.0, ... % seconds
        'starting_time_rate', params.downsample_fs, ... % Hz
        'data', lfp_all, ...
        'electrodes', electrode_table_region, ...
        'data_unit', 'uV');
    lfp = types.core.LFP('ElectricalSeries', lfp_electrical_series);
    ecephys_module = types.core.ProcessingModule(...
        'description', 'extracellular electrophysiology');
    ecephys_module.nwbdatainterface.set('LFP', lfp);
    nwb.processing.set('ecephys', ecephys_module);
    clear lfp_all

    % Estimate interconnectivity of single units
    if params.cnx

    end

    % Record eye trace
    if params.eye

    end

    % Save to NWB
    nwbExport(nwb, [out_file_path_tt 'sub-' recording_info.Subject '_ses-' recording_info.Session '.nwb']);
    disp(['SUCCESSFULLY SAVED: ' out_file_path_itt 'sub-' recording_info.Subject '_ses-' recording_info.Session '.nwb'])

    % Increment counter
    n_procd = n_procd + 1;

end

disp(['SUCCESSFULLY PROCESSED ' n_procd ' FILES.'])

end