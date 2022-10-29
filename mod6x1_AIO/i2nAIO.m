function nwb = i2nAIO(pp, nwb, recdev)

% in case of missing data...
default_adc_mapping = {'eye_x', 'eye_y', 'eye_pupil'};

downsample_factor = recdev.sampling_rate / 1000;

if strcmp(recdev.board_adc_channels(1).custom_channel_name, "ANALOG-IN-1")
    adc_map = convertCharsToStrings(default_adc_mapping);
    adc_map = adc_map(1:size(recdev.board_adc_data,1));
else
    for jj = 1 : numel(recdev.board_adc_channels)
        adc_map(jj) = convertCharsToStrings(recdev.board_adc_channels(jj).custom_channel_name);
    end
end

for jj = 1 : numel(adc_map)

    temp_dat = [];
    temp_dat(1,:) = downsample(recdev.board_adc_data(jj, :), downsample_factor);

    if strcmp(lower(adc_map(jj)), 'eye_x')

        find_y = find(ismember(lower(adc_map), "eye_y"));
        temp_dat(2,:) = downsample(recdev.board_adc_data(find_y, :), downsample_factor);

        find_p = find(ismember(lower(adc_map), "eye_pupil"));
        temp_pdat(1,:) = downsample(recdev.board_adc_data(find_p, :), downsample_factor);

        eye_position = types.core.SpatialSeries( ...
            'description', 'The position of the eye. Actual sampling rate = 500 Hz (Reported=1kHz)', ...
            'data', temp_dat, ...
            'starting_time_rate', 1000, ... % Hz
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        eye_tracking = types.core.EyeTracking();
        eye_tracking.spatialseries.set('eye_tracking', eye_position);

        pupil_diameter = types.core.TimeSeries( ...
            'description', 'Pupil diameter.', ...
            'data', temp_pdat, ...
            'starting_time_rate', 1000, ... % Hz
            'data_unit', 'arbitrary units', ...
            'timestamps', recdev.time_stamps_s_ds, ...
            'timestamps_unit', 'seconds' ...
            );

        pupil_tracking = types.core.PupilTracking();
        pupil_tracking.timeseries.set('pupil_diameter', pupil_diameter);

        nwb.acquisition.set('EyeTracking', eye_tracking);
        nwb.acquisition.set('PupilTracking', pupil_tracking);

        clear temp_* find_*
    end
end

nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

end