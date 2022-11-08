function updateNWBs(nwb_dirs)

add_convolutions = 0;
fix_event_codes = 0;

for ii = 1 : numel(nwb_dirs)
    dir_info = dir(nwb_dirs{ii});
    for jj = 1 : numel(dir_info)
        if strcmp(dir_info(jj).name(1), '.')
            continue
        end

        nwb = nwbRead([dir_info(jj).folder filesep dir_info(jj).name]);

        if fix_event_codes

            nwb.intervals = [];
            pp = pipelinePaths();
            recdev = readIntanHeader(paren(findDir('Z:\_DATA\_BL_DATA_PIPELINE\_0_RAW_DATA\', nwb.identifier),1));
            nwb = i2nDIO(pp, nwb, recdev);

        end

        if add_convolutions

            % ADD SPIKE TRAIN CONVOLUTIONS
            try nwb.processing.get('probe_0_suac').electricalseries.data
            catch
                try
                    conv_data = zeros(numel(nwb.units.spike_times_index.data(:)), size(nwb.acquisition.get('probe_0_lfp').electricalseries.get('probe_0_lfp_data').data(:,:), 2), 'single');
                catch
                    conv_data = zeros(numel(nwb.units.spike_times_index.data(:)), ceil(max(nwb.units.spike_times_index.data(:))*1000)+1000, 'single');
                end

                spike_times_indices = zeros(numel(nwb.units.spike_times.data(:)),1);
                for ii = 1 : numel(nwb.units.spike_times_index.data(:))
                    spike_times_indices(1:nwb.units.spike_times_index.data(ii)) = spike_times_indices(1:nwb.units.spike_times_index.data(ii)) + 1;
                end

                conv_data(sub2ind(size(conv_data), spike_times_indices, round(nwb.units.spike_times.data(:)*1000)))   = 1;

                Half_BW = ceil( 20 * 8 );
                x = 0 : Half_BW;
                k = [ zeros( 1, Half_BW ), ...
                    ( 1 - ( exp( -( x ./ 1 ) ) ) ) .* ( exp( -( x ./ 20) ) ) ];

                cnv_pre = mean(conv_data(:,1:floor(length(k)/2)),2)*ones(1,floor(length(k)/2));
                cnv_post = mean(conv_data(:,length(conv_data)-floor(length(k)/2):length(conv_data)),2)*ones(1,floor(length(k)/2));
                conv_data = conv2([ cnv_pre conv_data cnv_post ], k, 'valid') .* 1000;

                electrode_table_region_temp = types.hdmf_common.DynamicTableRegion( ...
                    'table', types.untyped.ObjectView(nwb.general_extracellular_ephys_electrodes), ...
                    'description', 'convolution peak channel references', ...
                    'data', nwb.units.vectordata.get('peak_channel_id').data(:));

                convolution_electrical_series = types.core.ElectricalSeries( ...
                    'electrodes', electrode_table_region_temp, ...
                    'starting_time', 0.0, ... % seconds
                    'starting_time_rate', 1000, ... % Hz
                    'data', conv_data, ...
                    'data_unit', 'spikes/second', ...
                    'filtering', 'Excitatory postsynaptic potential type convolution of spike rasters. kWidth=20', ...
                    'timestamps', nwb.acquisition.get('probe_0_lfp').electricalseries.get('probe_0_lfp_data').timestamps(:));

                suac_series = types.core.ProcessingModule('convolved_spike_train_data', convolution_electrical_series, ...
                    'description', 'Single units rasters convolved using EPSP kernel');
                nwb.processing.set('convolved_spike_train', suac_series);
            end

            nwbExport(nwb, [dir_info(jj).folder filesep dir_info(jj).name]);
            disp('done')
        end

    end
end

end