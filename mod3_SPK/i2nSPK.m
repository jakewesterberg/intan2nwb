function nwb = i2nSPK(pp, nwb, recdev, probe)

% ADD CODE TO CHECK FOR SPIKING ACTIVITY FROM DIFFERENT PROBES AND
% CONCATENATE. CURRENTLY THIS ONLY WORKS FOR SESSIONS WITH A SINGLE PROBE.
% ALSO WILL FAIL TO RUN IF A NWB FILE WITH UNIT DATA ALREADY EXISTS. SO WE
% START BY CLEARING THE UNITS FIELD

nwb.units = [];

reset(gpuDevice)

%Setup kilosort dirs
spk_file_path_itt = [pp.SPK_DATA nwb.identifier filesep 'probe-' ...
    num2str(probe.num) filesep]; % the raw data binary file is in this folder

if ~exist(spk_file_path_itt, 'dir')
    mkdir(spk_file_path_itt)
end

ops.chanMap = [pp.REPO 'forked_toolboxes' filesep 'Kilosort-2.0' ...
    filesep 'configFiles' filesep probe.type, '_kilosortChanMap.mat'];
run([pp.REPO 'forked_toolboxes' filesep 'Kilosort-2.0' ...
    filesep 'configFiles' filesep probe.type '_config.m'])

ops.trange      = [0 Inf]; % time range to sort
ops.NchanTOT    = probe.num_channels; % total number of channels in your recording

ops.fig = 0;
ops.fs = recdev.sampling_rate;

if ~exist(pp.SCRATCH, 'dir')
    mkdir(pp.SCRATCH)
end
ops.fproc = fullfile(pp.SCRATCH, 'temp_wh.dat'); % proc file on a fast SSD

% find the binary file
ops.fbinary = [pp.BIN_DATA nwb.identifier filesep nwb.identifier '_probe-' num2str(probe.num) '.bin'];

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(spk_file_path_itt, 'rez.mat'), 'rez', '-v7.3', '-nocompression');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, spk_file_path_itt);

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(spk_file_path_itt, 'rez2.mat');
save(fname, 'rez', '-v7.3', '-nocompression');

reset(gpuDevice)

% create json
json_struct = struct();

json_struct.directories.kilosort_output_directory = ...
    strrep(spk_file_path_itt, filesep, [filesep filesep]);

json_struct.waveform_metrics.waveform_metrics_file = ...
    strrep([spk_file_path_itt 'waveform_metrics.csv'], filesep, [filesep filesep]);

json_struct.ephys_params.sample_rate = recdev.sampling_rate;
json_struct.ephys_params.bit_volts = 0.195;
json_struct.ephys_params.num_channels = numel(rez.ycoords);
json_struct.ephys_params.reference_channels = []; %probe.num_channels/2;
json_struct.ephys_params.vertical_site_spacing = mean(diff(rez.ycoords));
json_struct.ephys_params.ap_band_file = ...
    strrep([pp.BIN_DATA nwb.identifier filesep nwb.identifier '_probe-' num2str(probe.num) '.bin'], filesep, [filesep filesep]);
json_struct.ephys_params.cluster_group_file_name = 'cluster_group.tsv.v2';
json_struct.ephys_params.reorder_lfp_channels = true;
json_struct.ephys_params.lfp_sample_rate = probe.downsample_fs;
json_struct.ephys_params.probe_type = probe.type;

json_struct.ks_postprocessing_params.within_unit_overlap_window = 0.000166;
json_struct.ks_postprocessing_params.between_unit_overlap_window = 0.000166;
json_struct.ks_postprocessing_params.between_unit_overlap_distance = 5;

json_struct.mean_waveform_params.mean_waveforms_file = ...
    strrep([spk_file_path_itt 'mean_waveforms.npy'], filesep, [filesep filesep]);
json_struct.mean_waveform_params.samples_per_spike = 82;
json_struct.mean_waveform_params.pre_samples = 20;
json_struct.mean_waveform_params.num_epochs = 1;
json_struct.mean_waveform_params.spikes_per_epoch = 1000;
json_struct.mean_waveform_params.spread_threshold = 0.12;
json_struct.mean_waveform_params.site_range = 16;

json_struct.noise_waveform_params.classifier_path = ...
    strrep([pp.REPO 'forked_toolboxes\ecephys_spike_sorting\modules\noise_templates\rf_classifier.pkl'], filesep, [filesep filesep]);
json_struct.noise_waveform_params.multiprocessing_worker_count = 10;

json_struct.quality_metrics_params.isi_threshold = 0.0015;
json_struct.quality_metrics_params.min_isi = 0.000166;
json_struct.quality_metrics_params.num_channels_to_compare = 7;
json_struct.quality_metrics_params.max_spikes_for_unit = 500;
json_struct.quality_metrics_params.max_spikes_for_nn = 10000;
json_struct.quality_metrics_params.n_neighbors = 4;
json_struct.quality_metrics_params.n_silhouette = 10000;
json_struct.quality_metrics_params.quality_metrics_output_file = ...
    strrep([spk_file_path_itt 'metrics_test.csv'], filesep, [filesep filesep]);
json_struct.quality_metrics_params.drift_metrics_interval_s = 51;
json_struct.quality_metrics_params.drift_metrics_min_spikes_per_interval = 10;
json_struct.quality_metrics_params.include_pc_metrics = true;

encodedJSON = jsonencode(json_struct);

fid = fopen([spk_file_path_itt 'ecephys_spike_sorting_input.json'], 'w');
fprintf(fid, encodedJSON);
fclose('all');

clear encodedJSON json_struct

fid = fopen([spk_file_path_itt 'ecephys_spike_sorting_adapter.bat'], 'w');
fprintf(fid, '%s\n', '@echo OFF');
fprintf(fid, '%s\n', ['set CONDAPATH=' pp.CONDA]);
fprintf(fid, '%s\n', 'set ENVNAME=ecephys');
fprintf(fid, '%s\n', 'if %ENVNAME%==base (set ENVPATH=%CONDAPATH%) else (set ENVPATH=%CONDAPATH%\envs\%ENVNAME%)');
fprintf(fid, '%s\n', 'call %CONDAPATH%\Scripts\activate.bat %ENVPATH%');
fprintf(fid, '%s\n', 'set GIT_PYTHON_REFRESH=quiet');
fprintf(fid, '%s\n', 'set PYTHONIOENCODING=utf-8');
fprintf(fid, '%s\n', ['cd ' pp.REPO 'forked_toolboxes\ecephys_spike_sorting']);
fprintf(fid, '%s\n', ['python -m ecephys_spike_sorting.modules.kilosort_postprocessing --input_json ' ...
    spk_file_path_itt 'ecephys_spike_sorting_input.json --output_json ' spk_file_path_itt 'ecephys_spike_sorting_kspp_output.json']);
fprintf(fid, '%s\n', ['python -m ecephys_spike_sorting.modules.mean_waveforms --input_json ' ...
    spk_file_path_itt 'ecephys_spike_sorting_input.json --output_json ' spk_file_path_itt 'ecephys_spike_sorting_waveforms_output.json']);
fprintf(fid, '%s\n', ['python -m ecephys_spike_sorting.modules.noise_templates --input_json ' ...
    spk_file_path_itt 'ecephys_spike_sorting_input.json --output_json ' spk_file_path_itt 'ecephys_spike_sorting_noise_output.json']);
fprintf(fid, '%s\n', ['python -m ecephys_spike_sorting.modules.quality_metrics --input_json '...
    spk_file_path_itt 'ecephys_spike_sorting_input.json --output_json ' spk_file_path_itt 'ecephys_spike_sorting_quality_output.json']);
fprintf(fid, '%s\n', 'call conda deactivate');
fclose('all');

system([spk_file_path_itt 'ecephys_spike_sorting_adapter.bat']);

unit_idents = unique(rez.st3(:,2))';
spike_times = cell(1, numel(unit_idents));
ctr_i=0;
for kk = unit_idents
    ctr_i = ctr_i + 1;
    spike_times{ctr_i} = (rez.st3(rez.st3(:,2)==kk,1)./recdev.sampling_rate).';

    % isi measures
    temp_isi = diff(spike_times{ctr_i});
    isi_mean(ctr_i) = mean(temp_isi);
    isi_cv(ctr_i) = std(temp_isi) / isi_mean(ctr_i);
    isi_0 = temp_isi(1:end-1);
    isi_1 = temp_isi(2:end);
    isi_lv(ctr_i) = (3/(numel(temp_isi)-1)) * sum(((isi_0-isi_1)./(isi_0+isi_1)).^2);

    clear isi_0 isi_1 temp_isi
end

% grab spike times and indices
[spike_times_vector, spike_times_index] = util.create_indexed_column(spike_times);

% grab the waveforms
mean_wave = readNPY([spk_file_path_itt 'mean_waveforms.npy']);
mean_wave_reshape = [];
for kk = 1 : numel(unit_idents)
    mean_wave_reshape = [mean_wave_reshape, squeeze(mean_wave(kk,:,:)).'];
end

clear mean_wave

% grab spike templates/amps
spike_amplitudes = readNPY([spk_file_path_itt 'amplitudes.npy']);
spike_amplitudes_index = spike_times_index.data(:);

% grab the metrics
fid = fopen([spk_file_path_itt 'metrics_test.csv'],'rt');
C = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', ',', 'HeaderLines', 1, 'EmptyValue', NaN);
fclose(fid);
[col_id, cluster_id, firing_rate, presence_ratio, isi_violations, amplitude_cutoff, isolation_distance, l_ratio, ...
    d_prime, nn_hit_rate, nn_miss_rate, silhouette_score, max_drift, cumulative_drift, epoch_name_quality_metrics, ...
    epoch_name_waveform_metrics, peak_channel_id, snr, waveform_duration, waveform_halfwidth, PT_ratio, repolarization_slope, ...
    recovery_slope, amplitude, spread, velocity_above, velocity_below] = deal(C{:});

clear C

% generate other indices
waveform_mean_index = probe.num_channels:probe.num_channels:probe.num_channels*numel(unit_idents);
local_index = 0:numel(unit_idents)-1; local_index = local_index';

% grab noise units
fid = fopen([spk_file_path_itt 'cluster_group.tsv.v2'],'rt');
C = textscan(fid, '%f %s', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
[quality_cluster_id, quality] = deal(C{:});
quality = quality(ismember(quality_cluster_id, cluster_id));

clear quality_cluster_id C

% gen colnames
colnames = {'snr';'cumulative_drift';'peak_channel_id';'quality';'local_index';'spread';'max_drift';'waveform_duration';'amplitude'; ...
    'amplitude_cutoff';'firing_rate';'nn_hit_rate';'nn_miss_rate';'silhouette_score';'isi_violations';'isolation_distance';'cluster_id'; ...
    'velocity_below';'repolarization_slope';'velocity_above';'l_ratio';'waveform_halfwidth';'presence_ratio';'PT_ratio';'recovery_slope';'d_prime'; ...
    'isi_mean';'isi_cv';'isi_lv';'waveform_mean_index';'spike_amplitudes_index'};

nwb.units = types.core.Units( ...
    'description', 'kilosorted and AllenSDK ecephys processed units', ...
    'colnames', colnames, ...
    'id', types.hdmf_common.ElementIdentifiers( ...
    'data', int64(0:numel(unit_idents) - 1)), ...
    'spike_times', spike_times_vector, ...
    'spike_times_index', spike_times_index, ...
    'waveform_mean', types.hdmf_common.VectorData('data', mean_wave_reshape, 'description', '2D Waveform'), ...
    'waveform_mean_index', types.hdmf_common.VectorData('data', waveform_mean_index, 'description', 'placeholder'), ...
    'PT_ratio', types.hdmf_common.VectorData('data', PT_ratio, 'description', 'placeholder'), ...
    'amplitude', types.hdmf_common.VectorData('data', amplitude, 'description', 'placeholder'), ...
    'amplitude_cutoff', types.hdmf_common.VectorData('data', amplitude_cutoff, 'description', 'placeholder'), ...
    'cluster_id', types.hdmf_common.VectorData('data', cluster_id, 'description', 'placeholder'), ...
    'cumulative_drift', types.hdmf_common.VectorData('data', cumulative_drift, 'description', 'placeholder'), ...
    'd_prime', types.hdmf_common.VectorData('data', d_prime, 'description', 'placeholder'), ...
    'firing_rate', types.hdmf_common.VectorData('data', firing_rate, 'description', 'placeholder'), ...
    'isi_violations', types.hdmf_common.VectorData('data', isi_violations, 'description', 'placeholder'), ...
    'isolation_distance', types.hdmf_common.VectorData('data', isolation_distance, 'description', 'placeholder'), ...
    'l_ratio', types.hdmf_common.VectorData('data', l_ratio, 'description', 'placeholder'), ...
    'local_index', types.hdmf_common.VectorData('data', local_index, 'description', 'placeholder'), ...
    'max_drift', types.hdmf_common.VectorData('data', max_drift, 'description', 'placeholder'), ...
    'nn_hit_rate', types.hdmf_common.VectorData('data', nn_hit_rate, 'description', 'placeholder'), ...
    'nn_miss_rate', types.hdmf_common.VectorData('data', nn_miss_rate, 'description', 'placeholder'), ...
    'peak_channel_id', types.hdmf_common.VectorData('data', peak_channel_id, 'description', 'placeholder'), ...
    'presence_ratio', types.hdmf_common.VectorData('data', presence_ratio, 'description', 'placeholder'), ...
    'quality', types.hdmf_common.VectorData('data', quality, 'description', 'placeholder'), ...
    'recovery_slope', types.hdmf_common.VectorData('data', recovery_slope, 'description', 'placeholder'), ...
    'repolarization_slope', types.hdmf_common.VectorData('data', repolarization_slope, 'description', 'placeholder'), ...
    'silhouette_score', types.hdmf_common.VectorData('data', silhouette_score, 'description', 'placeholder'), ...
    'snr', types.hdmf_common.VectorData('data', snr, 'description', 'placeholder'), ...
    'spike_amplitudes_index', types.hdmf_common.VectorData('data', spike_amplitudes_index, 'description', 'placeholder'), ...
    'spread', types.hdmf_common.VectorData('data', spread, 'description', 'placeholder'), ...
    'velocity_above', types.hdmf_common.VectorData('data', velocity_above, 'description', 'placeholder'), ...
    'velocity_below', types.hdmf_common.VectorData('data', velocity_below, 'description', 'placeholder'), ...
    'waveform_duration', types.hdmf_common.VectorData('data', waveform_duration, 'description', 'placeholder'), ...
    'waveform_halfwidth', types.hdmf_common.VectorData('data', waveform_halfwidth, 'description', 'placeholder'), ...
    'isi_mean', types.hdmf_common.VectorData('data', isi_mean, 'description', 'placeholder'), ...
    'isi_cv', types.hdmf_common.VectorData('data', isi_cv, 'description', 'placeholder'), ...
    'isi_lv', types.hdmf_common.VectorData('data', isi_lv, 'description', 'placeholder') ...
    );

clear isi_lv isi_cv isi_mean waveform_halfwidth waveform_duration velocity_below velocity_above spread spike_amplitdues_index
clear snr silhouette_score repolarization_slope recovery_slope quality presence_ratio peak_channel_id nn_miss_rate nn_hit_rate
clear max_drift local_index l_ratio isolation_distance isi_violations firing_rate d_prime cumulative_drift cluster_id amplitude_cutoff
clear amplitude PT_ratio waveform_mean_index mean_wave_reshape spike_times spike_times_index col_id epoch_name_quality_metrics
clear epoch_name_waveform_metrics

nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

nwb.units.vectordata.set('spike_amplitudes', types.hdmf_common.VectorData('description', 'placeholder', ...
    'data', spike_amplitudes));

clear spike_amplitudes

nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

%% CONVOLUTION
conv_data = zeros(numel(unit_idents), ceil(max(nwb.units.spike_times.data(:))*probe.downsample_fs)+probe.downsample_fs, 'single');

spike_times_indices = zeros(1, numel(nwb.units.spike_times.data(:)));
for ii = 1 : numel(unit_idents)
    spike_times_indices(1:nwb.units.spike_times_index.data(ii)) = spike_times_indices(1:nwb.units.spike_times_index.data(ii)) + 1;
end

conv_data(sub2ind(size(conv_data), spike_times_indices', round(nwb.units.spike_times.data(:)*probe.downsample_fs)))   = 1;

Half_BW = ceil( (20*(probe.downsample_fs/1000)) * 8 );
x = 0 : Half_BW;
k = [ zeros( 1, Half_BW ), ...
    ( 1 - ( exp( -( x ./ 1 ) ) ) ) .* ( exp( -( x ./ (probe.downsample_fs/1000)) ) ) ];

cnv_pre = mean(conv_data(:,1:floor(length(k)/2)),2)*ones(1,floor(length(k)/2));
cnv_post = mean(conv_data(:,length(conv_data)-floor(length(k)/2):length(conv_data)),2)*ones(1,floor(length(k)/2));
conv_data = conv2([ cnv_pre conv_data cnv_post ], k, 'valid') .* probe.downsample_fs;

electrode_table_region_temp = types.hdmf_common.DynamicTableRegion( ...
    'table', types.untyped.ObjectView(nwb.general_extracellular_ephys_electrodes), ...
    'description', 'convolution peak channel references', ...
    'data', nwb.units.vectordata.get('peak_channel_id').data(:));

convolution_electrical_series = types.core.ElectricalSeries( ...
    'electrodes', electrode_table_region_temp, ...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', probe.downsample_fs, ... % Hz
    'data', conv_data, ...
    'data_unit', 'spikes/second', ...
    'filtering', 'Excitatory postsynaptic potential type convolution of spike rasters. kWidth=20ms', ...
    'timestamps', (0:size(conv_data,2)-1)/probe.downsample_fs);

suac_series = types.core.ProcessingModule('convolved_spike_train_data', convolution_electrical_series, ...
    'description', 'Single units rasters convolved using EPSP kernel');
nwb.processing.set('convolved_spike_train', suac_series);

nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

end