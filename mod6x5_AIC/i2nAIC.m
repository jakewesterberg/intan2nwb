function nwb2 = i2nAIC(pp, nwb, recording_info, ii)

nwb2                                 = NwbFile;

% reformat existing
nwb2.identifier                      = recording_info.Identifier{ii};
nwb2.general_experimenter            = recording_info.Investigator{ii};
nwb2.general_session_id              = recording_info.Identifier{ii};
nwb2.general_experiment_description  = recording_info.Experiment_Description{ii};

nwb2.general_extracellular_ephys_electrodes = nwb.general_extracellular_ephys_electrodes;

% eye data
% eye_tracking = types.core.EyeTracking();
% pupil_tracking = types.core.PupilTracking();
% 
% eye_tracking.spatialseries.set('eye_1_tracking_data', nwb.acquisition.get('EyeTracking').eye_tracking);
% pupil_tracking.timeseries.set('pupil_1_diameter_data', nwb.acquisition.get('EyeTracking').pupil_tracking);
% 
% nwb.acquisition.set('eye_1_tracking', eye_tracking);
% nwb.acquisition.set('pupil_1_tracking', pupil_tracking);

% spiking additions
isi_mean                            = nan(1, numel(nwb.units.spike_times_index.data(:)));
isi_cv                              = nan(1, numel(nwb.units.spike_times_index.data(:)));
isi_lv                              = nan(1, numel(nwb.units.spike_times_index.data(:)));
stinds                              = [0; nwb.units.spike_times_index.data(:)];
unit_idents                         = 1:numel(nwb.units.spike_times_index.data(:));
ctr_i                               = 0;

for kk = unit_idents
    ctr_i = ctr_i + 1;

    % isi measures
    temp_isi = diff(nwb.units.spike_times.data(stinds(kk)+1:stinds(kk+1)));
    isi_mean(ctr_i) = mean(temp_isi);
    isi_cv(ctr_i) = std(temp_isi) / isi_mean(ctr_i);
    isi_0 = temp_isi(1:end-1);
    isi_1 = temp_isi(2:end);
    isi_lv(ctr_i) = (3/(numel(temp_isi)-1)) * sum(((isi_0-isi_1)./(isi_0+isi_1)).^2);

    clear isi_0 isi_1 temp_isi
end

% nwb2.units.vectordata.set('isi_mean', types.hdmf_common.VectorData('description', 'placeholder', 'data', isi_mean));
% nwb2.units.vectordata.set('isi_cv', types.hdmf_common.VectorData('description', 'placeholder', 'data', isi_cv));
% nwb2.units.vectordata.set('isi_lv', types.hdmf_common.VectorData('description', 'placeholder', 'data', isi_lv));

% single unit convolution
conv_data = zeros(numel(unit_idents), ceil(max(nwb.units.spike_times.data(:))*1000)+1000, 'single');

spike_times_indices = zeros(1, numel(nwb.units.spike_times.data(:)))-numel(unit_idents)-1;
for kk = 1 : numel(unit_idents)
    spike_times_indices(1:nwb.units.spike_times_index.data(kk)) = spike_times_indices(1:nwb.units.spike_times_index.data(kk)) + 1;
end

spike_times_indices = abs(spike_times_indices);
for kk = 1 : numel(unit_idents)
    conv_data(kk, round(nwb.units.spike_times.data(find(spike_times_indices==kk))*1000))   = 1;
end

rasters = int16(conv_data);

electrode_table_region_temp = types.hdmf_common.DynamicTableRegion( ...
    'table', types.untyped.ObjectView(nwb.general_extracellular_ephys_electrodes), ...
    'description', 'convolution peak channel references', ...
    'data', nwb.units.vectordata.get('peak_channel_id').data(:));

raster_electrical_series = types.core.ElectricalSeries( ...
    'electrodes', electrode_table_region_temp, ...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', 1000, ... % Hz
    'data', rasters, ...
    'data_unit', 'spikes', ...
    'filtering', 'spike times at discrete times', ...
    'timestamps', (0:size(conv_data,2)-1)/1000);

raster_series = types.core.ProcessingModule('spike_train_data', raster_electrical_series, ...
    'description', 'Spike trains in time');
nwb2.processing.set('spike_train', raster_series);

Half_BW = ceil( (20*(1000/1000)) * 8 );
x = 0 : Half_BW;
k = [ zeros( 1, Half_BW ), ...
    ( 1 - ( exp( -( x ./ 1 ) ) ) ) .* ( exp( -( x ./ (1000/1000)) ) ) ];
cnv_pre = mean(conv_data(:,1:floor(length(k)/2)),2)*ones(1,floor(length(k)/2));
cnv_post = mean(conv_data(:,length(conv_data)-floor(length(k)/2):length(conv_data)),2)*ones(1,floor(length(k)/2));

for mm = 1 : size(conv_data,1)
    conv_data(mm,:) = conv([cnv_pre(mm,:) conv_data(mm,:) cnv_post(mm,:)], k, 'valid') .* 1000;
end

convolution_electrical_series = types.core.ElectricalSeries( ...
    'electrodes', electrode_table_region_temp, ...
    'starting_time', 0.0, ... % seconds
    'starting_time_rate', 1000, ... % Hz
    'data', conv_data, ...
    'data_unit', 'spikes/second', ...
    'filtering', 'Excitatory postsynaptic potential type convolution of spike rasters. kWidth=20ms', ...
    'timestamps', (0:size(conv_data,2)-1)/1000);

suac_series = types.core.ProcessingModule('convolved_spike_train_data', convolution_electrical_series, ...
    'description', 'Single units rasters convolved using EPSP kernel');
nwb2.processing.set('convolved_spike_train', suac_series);

% add lfp data
raw_data_dir = findDir(pp.RAW_DATA, recording_info.Identifier{ii});
[~, dir_name_temp] = fileparts(raw_data_dir);
probe_files = findFiles([pp.RAW_DATA dir_name_temp filesep], 'probe');

p_ctr = 1;
for kk = 1 : numel(probe_files)

    no_luck = 1;
    while no_luck
        try
            nwb2.general_extracellular_ephys.set(['probe' alphabet(kk)], nwb.general_extracellular_ephys.get(['probe' alphabet(p_ctr)]));
            p_ctr = p_ctr + 1;
            no_luck = 0;
        catch
            p_ctr = p_ctr + 1;
        end
    end

    nwb_lfp = nwbRead(probe_files{kk});
    lfp_electrical_series = nwb_lfp.acquisition.get(['probe_' num2str(kk-1) '_lfp']).electricalseries.get(['probe_' num2str(kk-1) '_lfp_data']);
    lfp_electrical_series.timestamps = lfp_electrical_series.timestamps(:);

    lfp_electrical_series.data = lfp_electrical_series.data(:,:);

    lfp_series = types.core.LFP(['probe_' num2str(kk-1) '_lfp_data'], lfp_electrical_series);
    nwb2.acquisition.set(['probe_' num2str(kk-1) '_lfp'], lfp_series);
end

% event coding
if strcmp(nwb.general_stimulus, 'OpenScopeGlobalLocalOddball')
    event_data{1} = ALLENINSTITUTE_PassiveGLOv1(nwb);
    event_data{2} = ALLENINSTITUTE_RFMappingv1(nwb);
    event_data{3} = ALLENINSTITUTE_Optotaggingv1(nwb);
end

for jj = 1 : numel(event_data)
    temp_fields = fields(event_data{jj});
    temp_fields = temp_fields(~strcmp(temp_fields, 'task'));

    eval_str = [];
    for kk = 1 : numel(temp_fields)
        eval_str = ...
            [ eval_str ...
            ',convertStringsToChars("' ...
            temp_fields{kk} ...
            '"), types.hdmf_common.VectorData(convertStringsToChars("data"), event_data{jj}.' ...
            temp_fields{kk} ', convertStringsToChars("description"), convertStringsToChars("placeholder"))'];
    end
    eval_str = [
        'trials=types.core.TimeIntervals(convertStringsToChars("description"), convertStringsToChars("events"), convertStringsToChars("colnames"),temp_fields' ...
        eval_str ');']; ...

    eval(eval_str); clear eval_str
    nwb2.intervals.set(event_data{jj}.task, trials); clear trials
end

nwbExport(nwb2, [pp.NWB_DATA nwb2.identifier '.nwb']);

end