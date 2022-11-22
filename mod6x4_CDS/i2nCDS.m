function nwb = i2nCDS(pp, nwb, recdev, probe)

% Let's change this to working on the bin files?

% Initialize DC offset filter
[DC_offset_bwb, DC_offset_bwa] = butter(1, 0.1/(recdev.sampling_rate/2), 'high');

% Initialize filter information
[muae_bwb, muae_bwa] = butter(2, [500 5000]/(recdev.sampling_rate/2), 'bandpass');
[muae_power_bwb, muae_power_bwa] = butter(4, 250/(recdev.sampling_rate/2), 'low');
[lfp_bwb, lfp_bwa] = butter(2, [0.1 500]/(recdev.sampling_rate/2), 'bandpass');

% Load the correct channel map file
load([probe.type '.mat'], 'channel_map')

% Initialize data matrices. Need to fix for multiprobe
lfp = zeros(probe.num_channels, recdev.downsample_size);
muae = zeros(probe.num_channels, recdev.downsample_size);

% Set computations to CPU, you are limited by RAM/VRAM at this
% point. might as well use whatever you have more of...
test_fid = fopen(recdev.in_file_path + "\amp-" + recdev.amplifier_channels(1).native_channel_name + ".dat");
test_size = byteSize(double(fread(test_fid, probe.num_samples, 'int16')) .* 0.195);
%workers = floor((gpuDevice().AvailableMemory) / (6*test_size));
mem = memory;
workers = floor((mem.MemAvailableAllArrays) / (6*test_size)); % change six to higher number if mem issues
if workers > feature('numcores')
    workers = feature('numcores');
elseif workers == 0
    workers = 1;
end
fclose(test_fid);
clear test_size

pvar_amp_ch = cat(1,{recdev.amplifier_channels.native_channel_name});
pvar_ds_factor = probe.downsample_factor;

in_file_path = recdev.in_file_path;
num_samples = probe.num_samples;

if ~isempty(gcp('nocreate'))
    delete(gcp);
end
pool1 = parpool(workers);
parfor kk = 1:probe.num_channels
    % Open file and init data
    current_fid             = fopen(in_file_path + "\amp-" + pvar_amp_ch{kk} + ".dat");

    % Setup array on GPU or in mem depending on run parameters
    %current_data            = gpuArray(double(fread(current_fid, probe.num_samples, 'uint16')  - 32768) * 0.195);
    current_data            = double(fread(current_fid, num_samples, 'int16')) .* 0.195;

    % Do data type specific filtering
    %             muae(kk,:)  = gather(downsample(filtfilt(muae_power_bwb, muae_power_bwa, ...
    %                 abs(filtfilt(muae_bwb, muae_bwa, ...
    %                 filtfilt(DC_offset_bwb, DC_offset_bwa, ...
    %                 current_data)))), pvar_ds_factor));
    %             lfp(kk,:)   = gather(downsample(filtfilt(lfp_bwb, lfp_bwa, ...
    %                 filtfilt(DC_offset_bwb, DC_offset_bwa, ...
    %                 current_data)), pvar_ds_factor));
    %             reset(gpuDevice)
    muae(kk,:)  = downsample(filtfilt(muae_power_bwb, muae_power_bwa, ...
        abs(filtfilt(muae_bwb, muae_bwa, ...tha
        filtfilt(DC_offset_bwb, DC_offset_bwa, ...
        current_data)))), pvar_ds_factor);
    lfp(kk,:)   = downsample(filtfilt(lfp_bwb, lfp_bwa, ...
        current_data), pvar_ds_factor);
    %            reset(gpuDevice)

    % Close file
    fclose(current_fid);
    disp([num2str(kk) '/' num2str(probe.num_channels) ' COMPLETED.'])

end
delete(pool1)
clear pvar_*

%Rearrange the channels to the order on the probe (starts at 0, +1 so it
%matches matlab indexing)
muae = muae(channel_map+1,:);
lfp = lfp(channel_map+1,:);

lfp_electrical_series = types.core.ElectricalSeries( ...
    'electrodes', probe.electrode_table_region,...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', probe.downsample_fs, ... % Hz
    'data', lfp, ...
    'data_unit', 'uV', ...
    'filtering', '4th order Butterworth 1-250 Hz', ...
    'timestamps', recdev.time_stamps_s_ds);

lfp_series = types.core.LFP(['probe_' num2str(probe.num) '_lfp_data'], lfp_electrical_series);
nwb.acquisition.set(['probe_' num2str(probe.num) '_lfp'], lfp_series);
clear lfp*

muae_electrical_series = types.core.ElectricalSeries( ...
    'electrodes', probe.electrode_table_region,...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', probe.downsample_fs, ... % Hz
    'data', muae, ...
    'data_unit', 'uV', ...
    'filtering', '4th order Butterworth 500-500 Hz, full-wave rectified, then low pass 4th order Butterworth 250 Hz (DC offset high-pass 1st order Butterworth 0.1 Hz)', ...
    'timestamps', recdev.time_stamps_s_ds);

muae_series = types.core.LFP(['probe_' num2str(probe.num) '_muae_data'], muae_electrical_series);
nwb.acquisition.set(['probe_' num2str(probe.num) '_muae'], muae_series);
clear muae*

nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

end