function nwb = i2nAIO(nwb, recdev)

downsample_factor = recdev.sampling_rate / 1000;

adc_map = recdev.adc_map;

pd_ctr = 0;
rew_ctr = 0;
laser_ctr = 0;
for jj = 1 : numel(adc_map)

    temp_dat = [];
    temp_dat(1,:) = downsample(recdev.board_adc_data(jj, :), downsample_factor);

    if strcmp(lower(adc_map(jj)), 'eye_1_x')

        find_y = find(ismember(lower(adc_map), "eye_1_y"));
        temp_dat(2,:) = downsample(recdev.board_adc_data(find_y, :), downsample_factor);

        find_p = find(ismember(lower(adc_map), "eye_1_pupil"));
        temp_pdat(1,:) = downsample(recdev.board_adc_data(find_p, :), downsample_factor);

        eye_position = types.core.SpatialSeries( ...
            'description', 'The position of the eye. Actual sampling rate = 500 Hz (Reported=1kHz)', ...
            'data', temp_dat, ...
            'starting_time_rate', 1000, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        eye_tracking = types.core.EyeTracking();
        eye_tracking.spatialseries.set('eye_1_tracking_data', eye_position);

        pupil_diameter = types.core.TimeSeries( ...
            'description', 'Pupil diameter.', ...
            'data', temp_pdat, ...
            'starting_time_rate', 1000, ... % Hz
            'data_unit', 'arbitrary units', ...
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        pupil_tracking = types.core.PupilTracking();
        pupil_tracking.timeseries.set('pupil_1_diameter_data', pupil_diameter);

        nwb.acquisition.set('eye_1_tracking', eye_tracking);
        nwb.acquisition.set('pupil_1_tracking', pupil_tracking);

        clear temp_* find_*
    end

    if strcmp(lower(adc_map(jj)), 'eye_2_x')

        find_y = find(ismember(lower(adc_map), "eye_2_y"));
        temp_dat(2,:) = downsample(recdev.board_adc_data(find_y, :), downsample_factor);

        find_p = find(ismember(lower(adc_map), "eye_2_pupil"));
        temp_pdat(1,:) = downsample(recdev.board_adc_data(find_p, :), downsample_factor);

        eye_position = types.core.SpatialSeries( ...
            'description', 'The position of the eye. Actual sampling rate = 500 Hz (Reported=1kHz)', ...
            'data', temp_dat, ...
            'starting_time_rate', 1000, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        eye_tracking = types.core.EyeTracking();
        eye_tracking.spatialseries.set('eye_2_tracking_data', eye_position);

        pupil_diameter = types.core.TimeSeries( ...
            'description', 'Pupil diameter.', ...
            'data', temp_pdat, ...
            'starting_time_rate', 1000, ... % Hz
            'data_unit', 'arbitrary units', ...
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        pupil_tracking = types.core.PupilTracking();
        pupil_tracking.timeseries.set('pupil_2_diameter_data', pupil_diameter);

        nwb.acquisition.set('eye_2_tracking', eye_tracking);
        nwb.acquisition.set('pupil_2_tracking', pupil_tracking);

        clear temp_* find_*
    end

    if strcmp(lower(adc_map(jj)), 'photodiode')
        pd_ctr = pd_ctr + 1;

        pd_state = types.core.TimeSeries( ...
            'description', 'photodiode electrical state in V', ...
            'data', temp_dat, ...
            'data_unit', 'Volts', ...
            'starting_time_rate', 1000, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        pd_tracking = types.core.BehavioralTimeSeries();
        pd_tracking.timeseries.set(['photodiode_' num2str(pd_ctr) '_tracking_data'], pd_state);

        nwb.acquisition.set(['photodiode_' num2str(pd_ctr) '_tracking'], pd_tracking);

        [bwb, bwa] = butter(1, 30/1000, 'low');
        temp_data_2 = find(ischange(filtfilt(bwb, bwa, abs(temp_dat)))) / 1000;

        pd_dat_1 = types.hdmf_common.VectorData('data', temp_data_2, 'description', 'time in s'); clear temp_data_2

        trials = types.core.TimeIntervals('description', 'photodiode detected changes in seconds using ischange on mean', ...
            'start_time', pd_dat_1, ...
            'stop_time', pd_dat_1, ...
            'colnames', {'start_time', 'stop_time'});

        nwb.intervals.set(['photodiode_' num2str(pd_ctr) '_detected_changes'], trials); clear trials

    end

    if strcmp(lower(adc_map(jj)), 'reward')
        rew_ctr = rew_ctr + 1;

        rew_state = types.core.TimeSeries( ...
            'description', 'reward signal electrical state in V', ...
            'data', temp_dat, ...
            'data_unit', 'Volts', ...
            'starting_time_rate', 1000, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        rew_tracking = types.core.BehavioralTimeSeries();
        rew_tracking.timeseries.set(['reward_' num2str(rew_ctr) '_tracking_data'], rew_state);

        nwb.acquisition.set(['reward_' num2str(rew_ctr) '_tracking'], rew_tracking);

        temp_data_2 = find(ischange(abs(temp_dat))) / 1000;
        if mod(numel(temp_data_2),2) ~= 0
            warning('reward detection may be inaccurate')
            threshold_val = 1;
            inacc_size = size(temp_data_2, 2);
            while size(temp_data_2) ~= inacc_size-1
                if size(temp_data_2) <= inacc_size-2
                    threshold_val = threshold_val-0.0001;
                    temp_data_2 = find(ischange(abs(temp_dat), 'Threshold', threshold_val)) / 1000;
                else
                    threshold_val = threshold_val+0.0005;
                    temp_data_2 = find(ischange(abs(temp_dat), 'Threshold', threshold_val)) / 1000;
                end
            end
        end
        temp_data_3 = temp_data_2(1:2:end);
        temp_data_4 = temp_data_2(2:2:end); clear temp_data_2

        rew_dat_1 = types.hdmf_common.VectorData('data', temp_data_3, 'description', 'time in s'); clear temp_data_3
        rew_dat_2 = types.hdmf_common.VectorData('data', temp_data_4, 'description', 'time in s'); clear temp_data_4

        trials = types.core.TimeIntervals('description', 'reward detected in seconds using ischange on mean', ...
            'start_time', rew_dat_1, ...
            'stop_time', rew_dat_2, ...
            'colnames', {'start_time', 'stop_time'});

        nwb.intervals.set(['reward_' num2str(rew_ctr) '_detected'], trials); clear trials

    end

    if strcmp(lower(adc_map(jj)), 'laser')
        laser_ctr = laser_ctr + 1;

        laser_state = types.core.TimeSeries( ...
            'description', 'reward signal electrical state in V', ...
            'data', temp_dat, ...
            'data_unit', 'Volts', ...
            'starting_time_rate', 1000, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        laser_tracking = types.core.BehavioralTimeSeries();
        laser_tracking.timeseries.set(['laser_' num2str(laser_ctr) '_tracking_data'], laser_state);

        nwb.acquisition.set(['laser_' num2str(laser_ctr) '_tracking'], laser_tracking);

        % SPECIFIC TO THE BASTOS LAB LASER SETUP
        laser_on = temp_dat >= 0.7;
        laser_off = temp_dat < 0.7;

        laser_peaks = intersect( ...
            find(islocalmax(temp_dat, 'MinProminence', 0.1)), ...
            find(laser_on));

        laser_peaks_volts = round(temp_dat(laser_peaks), 1);
        laser_peaks_bin = kmeansElbow([laser_peaks_volts; laser_peaks_volts]');
        K = unique(laser_peaks_bin);
        for ls = K'
            laser_peaks_volts(laser_peaks_bin == ls) = ...
                round(mean(laser_peaks_volts(laser_peaks_bin == ls)), 1);
        end
        laser_peaks_volts = double(laser_peaks_volts( ...
            [1, (find(diff(laser_peaks) > 500) + 1)]))';

        end_laser_stim = laser_peaks(find(diff(laser_peaks) > 500));
        for ss = 1 : numel(end_laser_stim)
            end_laser_stim(ss) = end_laser_stim(ss) + ...
                find(laser_off(end_laser_stim(ss):end)==1, 1, "first");
        end

        start_laser_stim = laser_peaks(find(diff(laser_peaks) > 500) + 1);
        for ss = 1 : numel(start_laser_stim)
            start_laser_stim(ss) = ...
                find(laser_off(1:start_laser_stim(ss))==1, 1, "last");
        end

        end_laser_stim = [end_laser_stim find(laser_on==1, 1, "last")]';
        start_laser_stim = [find(laser_on==1, 1, "first") start_laser_stim]';

        frequency_ct = nan(numel(start_laser_stim), 1);
        for ls = 1 : numel(start_laser_stim)
            frequency_ct(ls) = sum( ...
                laser_peaks < end_laser_stim(ls) & ...
                laser_peaks > start_laser_stim(ls));
        end
        frequency_bin = kmeansElbow( ...
            [frequency_ct frequency_ct]);

        K = unique(frequency_bin);
        for ls = K'
            frequency_ct(frequency_bin == ls) = ...
                round(mean(frequency_ct(frequency_bin == ls)));
        end

        start_laser_stim = start_laser_stim ./ 1000;
        end_laser_stim = end_laser_stim ./ 1000;

        laser_type_id = {};
        for ls = 1 : numel(frequency_ct)
            laser_type_id{ls,1} = [num2str(frequency_ct(ls)) '_hz_pulse_train'];
        end

        laser_start = types.hdmf_common.VectorData('data', start_laser_stim, 'description', 'time in s');
        laser_stop = types.hdmf_common.VectorData('data', end_laser_stim, 'description', 'time in s');

        laser_level = types.hdmf_common.VectorData('data', double(laser_peaks_volts)', 'description', 'stim au level');
        laser_frequency = types.hdmf_common.VectorData('data', frequency_ct, 'description', 'stim sine freq');
        laser_type = types.hdmf_common.VectorData('data', laser_type_id, 'description', 'stim label');

        trials = types.core.TimeIntervals('description', 'detected laser events from analog input', ...
            'start_time', laser_start, ...
            'stop_time', laser_stop, ...
            'level', laser_level, ...
            'frequency', laser_frequency, ...
            'type', laser_type, ...
            'colnames', {'start_time', 'stop_time', 'level', 'frequency'});

        nwb.intervals.set(['laser_' num2str(rew_ctr) '_empirical_decode'], trials); clear trials

    end
end

%% Add saccade detection once eye data is established.
% if any(strcmp(lower(adc_map(:)), 'eye_1_x'))
%
%     otero_input = [ ...
%         nwb.acquisition.get('eye_1_tracking').spatialseries.get('eye_1_tracking_data').timestamps(:), ...
%         nwb.acquisition.get('eye_1_tracking').spatialseries.get('eye_1_tracking_data').data(1,:)', ...
%         nwb.acquisition.get('eye_1_tracking').spatialseries.get('eye_1_tracking_data').data(2,:)'];
%
%     if any(strcmp(lower(adc_map(:)), 'eye_2_x'))
%         otero_input = [otero_input, ...
%             nwb.acquisition.get('eye_2_tracking').spatialseries.get('eye_2_tracking_data').data(1,:)', ...
%             nwb.acquisition.get('eye_2_tracking').spatialseries.get('eye_2_tracking_data').data(2,:)'];
%     else
%         otero_input = [otero_input, ...
%             nan(size(otero_input,1), 1), ...
%             nan(size(otero_input,1), 1)];
%     end
%
%     % blink detection is not implemented. note that blinks might be
%     % saccades in instances where a trial was aborted.
%
%     recording           = ClusterDetection.EyeMovRecording.Create( ...
%         pp.SCRATCH, ...
%         'temp_sac', ...
%         otero_input, ...
%         zeros(size(otero_input,1), 1), ...
%         1000);
%
%     [saccades, stats]    = recording.FindSaccades();
% end
%nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

end