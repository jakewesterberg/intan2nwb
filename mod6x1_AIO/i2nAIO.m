function nwb = i2nAIO(pp, nwb, recdev)

downsample_factor = recdev.sampling_rate / 1250;

adc_map = recdev.adc_map;

pd_ctr = 0;
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
            'starting_time_rate', 1250, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        eye_tracking = types.core.EyeTracking();
        eye_tracking.spatialseries.set('eye_1_tracking_data', eye_position);

        pupil_diameter = types.core.TimeSeries( ...
            'description', 'Pupil diameter.', ...
            'data', temp_pdat, ...
            'starting_time_rate', 1250, ... % Hz
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
            'starting_time_rate', 1250, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        eye_tracking = types.core.EyeTracking();
        eye_tracking.spatialseries.set('eye_2_tracking_data', eye_position);

        pupil_diameter = types.core.TimeSeries( ...
            'description', 'Pupil diameter.', ...
            'data', temp_pdat, ...
            'starting_time_rate', 1250, ... % Hz
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
            'starting_time_rate', 1250, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        pd_tracking = types.core.BehavioralTimeSeries();
        pd_tracking.timeseries.set(['photodiode_' num2str(pd_ctr) '_tracking_data'], pd_state);

        nwb.acquisition.set(['photodiode_' num2str(pd_ctr) '_tracking'], pd_tracking);

        [bwb, bwa] = butter(1, 30/1250, 'low');
        temp_data_2 = find(ischange(filtfilt(bwb, bwa, abs(temp_dat)))) / 1250;

        pd_dat_1 = types.hdmf_common.VectorData('data', temp_data_2, 'description', 'time in s'); clear temp_data_2

        trials = types.core.TimeIntervals('description', 'photodiode detected changes in seconds using ischange on mean', ...
            'start_time', pd_dat_1, ...
            'stop_time', pd_dat_1, ...
            'colnames', {'start_time', 'stop_time'});

        nwb.intervals.set(['photodiode_' num2str(pd_ctr) '_detected_changes'], trials); clear trials

    end

    if strcmp(lower(adc_map(jj)), 'reward')

        rew_state = types.core.TimeSeries( ...
            'description', 'reward signal electrical state in V', ...
            'data', temp_dat, ...
            'data_unit', 'Volts', ...
            'starting_time_rate', 1250, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        rew_tracking = types.core.BehavioralTimeSeries();
        rew_tracking.timeseries.set(['reward_' num2str(pd_ctr) '_tracking_data'], rew_state);

        nwb.acquisition.set(['reward_' num2str(pd_ctr) '_tracking'], rew_tracking);

        temp_data_2 = find(ischange(abs(temp_dat))) / 1250;
        if mod(numel(temp_data_2),2) ~= 0
            warning('reward detection may be inaccurate')
        end
        temp_data_3 = temp_data_2(1:2:end);
        temp_data_4 = temp_data_2(2:2:end); clear temp_data_2

        rew_dat_1 = types.hdmf_common.VectorData('data', temp_data_3, 'description', 'time in s'); clear temp_data_3
        rew_dat_2 = types.hdmf_common.VectorData('data', temp_data_4, 'description', 'time in s'); clear temp_data_4

        trials = types.core.TimeIntervals('description', 'reward detected in seconds using ischange on mean', ...
            'start_time', rew_dat_1, ...
            'stop_time', rew_dat_2, ...
            'colnames', {'start_time', 'stop_time'});

        nwb.intervals.set(['reward_' num2str(pd_ctr) '_detected'], trials); clear trials

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
%         1250);
% 
%     [saccades, stats]    = recording.FindSaccades();
% end
%nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

end