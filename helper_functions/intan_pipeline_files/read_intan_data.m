%Function that when used, you can either put the path to the folder that
%contains your probe recording data (has to be the "one file per channel"
%format from the intan recording system), or you can just run it, and
%select the folder yourself from the file explorer on your system, 
function probe_data = read_intan_data(folder_path,desired_port)
    %If a valid path was not provided in the input, this opens up a UI for
    %the user to select the folder that contains the probe data themselves.
    if ~exist('folder_path','var')
        folder_path = uigetdir;
        if isnumeric(folder_path)
            disp("Folder not selected")
            return;
        end
    end

    %If the user wishes to analyze one port at a time (considering if we
    %did multiple probes at once, it might be impractical to handle all of
    %the channels at once)
    %Then input the port value (Usually A-H)
    %However, if all of them are desired, enter 'ALL', or leave it blank,
    %and all of the channels will be read into memory
    if ~exist('desired_port','var')
        desired_port = input("Which probe to read? (A - H (ALL)): ","s");
        desired_port = char(desired_port);
        desired_port = upper(desired_port);
    end
    
    %This rhd file contains information about the probe
    info_fid = fopen(folder_path + "\info.rhd");
    
    %Unable to find an info file
    if info_fid == -1
        if (fopen(folder_path + "settings.xml"))
            %Just to aid in debugging, its possible to record in a
            %different data format than we were planning to use, this
            %message provides a better hint.
            disp("Make sure to use 'One File Per Channel', this function cannot read the other file formats")
        else
            disp("Could not locate info.rhd file, check if you have selected the folder containing your probe recording")
        end
        return
    end

    %% Taken directly from read_Intan_RHD2000_file.m, until this section ends

    %This is a ready made solution to read the cryptic data bits in
    %"info.rhd", which describes important stuff like, number of channels,
    %sampling rate, etc etc 

    %The only thing I changed was naming it "info_fid" so I know which file
    %it is that I'm opening
    
    % Check 'magic number' at beginning of file to make sure this is an Intan
    % Technologies RHD2000 data file.
    magic_number = fread(info_fid, 1, 'uint32');
    if magic_number ~= hex2dec('c6912702')
        error('Unrecognized file type.');
    end
    
    % Read version number.
    data_file_main_version_number = fread(info_fid, 1, 'int16');
    data_file_secondary_version_number = fread(info_fid, 1, 'int16');
    
    fprintf(1, '\n');
    fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
        data_file_main_version_number, data_file_secondary_version_number);
    fprintf(1, '\n');
    
    if (data_file_main_version_number == 1)
        num_samples_per_data_block = 60;
    else
        num_samples_per_data_block = 128;
    end
    
    % Read information of sampling rate and amplifier frequency settings.
    sample_rate = fread(info_fid, 1, 'single');
    dsp_enabled = fread(info_fid, 1, 'int16');
    actual_dsp_cutoff_frequency = fread(info_fid, 1, 'single');
    actual_lower_bandwidth = fread(info_fid, 1, 'single');
    actual_upper_bandwidth = fread(info_fid, 1, 'single');
    
    desired_dsp_cutoff_frequency = fread(info_fid, 1, 'single');
    desired_lower_bandwidth = fread(info_fid, 1, 'single');
    desired_upper_bandwidth = fread(info_fid, 1, 'single');
    
    % This tells us if a software 50/60 Hz notch filter was enabled during
    % the data acquisition.
    notch_filter_mode = fread(info_fid, 1, 'int16');
    notch_filter_frequency = 0;
    if (notch_filter_mode == 1)
        notch_filter_frequency = 50;
    elseif (notch_filter_mode == 2)
        notch_filter_frequency = 60;
    end
    
    desired_impedance_test_frequency = fread(info_fid, 1, 'single');
    actual_impedance_test_frequency = fread(info_fid, 1, 'single');
    
    % Place notes in data strucure
    notes = struct( ...
        'note1', fread_QString(info_fid), ...
        'note2', fread_QString(info_fid), ...
        'note3', fread_QString(info_fid) );
        
    % If data file is from GUI v1.1 or later, see if temperature sensor data
    % was saved.
    num_temp_sensor_channels = 0;
    if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 1) ...
        || (data_file_main_version_number > 1))
        num_temp_sensor_channels = fread(info_fid, 1, 'int16');
    end
    
    % If data file is from GUI v1.3 or later, load board mode.
    board_mode = 0;
    if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 3) ...
        || (data_file_main_version_number > 1))
        board_mode = fread(info_fid, 1, 'int16');
    end
    
    % If data file is from v2.0 or later (Intan Recording Controller),
    % load name of digital reference channel.
    if (data_file_main_version_number > 1)
        reference_channel = fread_QString(info_fid);
    end
    
    % Place frequency-related information in data structure.
    frequency_parameters = struct( ...
        'amplifier_sample_rate', sample_rate, ...
        'aux_input_sample_rate', sample_rate / 4, ...
        'supply_voltage_sample_rate', sample_rate / num_samples_per_data_block, ...
        'board_adc_sample_rate', sample_rate, ...
        'board_dig_in_sample_rate', sample_rate, ...
        'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
        'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
        'dsp_enabled', dsp_enabled, ...
        'desired_lower_bandwidth', desired_lower_bandwidth, ...
        'actual_lower_bandwidth', actual_lower_bandwidth, ...
        'desired_upper_bandwidth', desired_upper_bandwidth, ...
        'actual_upper_bandwidth', actual_upper_bandwidth, ...
        'notch_filter_frequency', notch_filter_frequency, ...
        'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
        'actual_impedance_test_frequency', actual_impedance_test_frequency );
    
    % Define data structure for spike trigger settings.
    spike_trigger_struct = struct( ...
        'voltage_trigger_mode', {}, ...
        'voltage_threshold', {}, ...
        'digital_trigger_channel', {}, ...
        'digital_edge_polarity', {} );
    
    new_trigger_channel = struct(spike_trigger_struct);
    spike_triggers = struct(spike_trigger_struct);
    
    % Define data structure for data channels.
    channel_struct = struct( ...
        'native_channel_name', {}, ...
        'custom_channel_name', {}, ...
        'native_order', {}, ...
        'custom_order', {}, ...
        'board_stream', {}, ...
        'chip_channel', {}, ...
        'port_name', {}, ...
        'port_prefix', {}, ...
        'port_number', {}, ...
        'electrode_impedance_magnitude', {}, ...
        'electrode_impedance_phase', {} );
    
    new_channel = struct(channel_struct);
    
    % Create structure arrays for each type of data channel.
    amplifier_channels = struct(channel_struct);
    aux_input_channels = struct(channel_struct);
    supply_voltage_channels = struct(channel_struct);
    board_adc_channels = struct(channel_struct);
    board_dig_in_channels = struct(channel_struct);
    board_dig_out_channels = struct(channel_struct);
    
    amplifier_index = 1;
    aux_input_index = 1;
    supply_voltage_index = 1;
    board_adc_index = 1;
    board_dig_in_index = 1;
    board_dig_out_index = 1;
    
    % Read signal summary from data file header.
    
    number_of_signal_groups = fread(info_fid, 1, 'int16');
    
    for signal_group = 1:number_of_signal_groups
        signal_group_name = fread_QString(info_fid);
        signal_group_prefix = fread_QString(info_fid);
        signal_group_enabled = fread(info_fid, 1, 'int16');
        signal_group_num_channels = fread(info_fid, 1, 'int16');
        signal_group_num_amp_channels = fread(info_fid, 1, 'int16');
    
        if (signal_group_num_channels > 0 && signal_group_enabled > 0)
            new_channel(1).port_name = signal_group_name;
            new_channel(1).port_prefix = signal_group_prefix;
            new_channel(1).port_number = signal_group;
            for signal_channel = 1:signal_group_num_channels
                new_channel(1).native_channel_name = fread_QString(info_fid);
                new_channel(1).custom_channel_name = fread_QString(info_fid);
                new_channel(1).native_order = fread(info_fid, 1, 'int16');
                new_channel(1).custom_order = fread(info_fid, 1, 'int16');
                signal_type = fread(info_fid, 1, 'int16');
                channel_enabled = fread(info_fid, 1, 'int16');
                new_channel(1).chip_channel = fread(info_fid, 1, 'int16');
                new_channel(1).board_stream = fread(info_fid, 1, 'int16');
                new_trigger_channel(1).voltage_trigger_mode = fread(info_fid, 1, 'int16');
                new_trigger_channel(1).voltage_threshold = fread(info_fid, 1, 'int16');
                new_trigger_channel(1).digital_trigger_channel = fread(info_fid, 1, 'int16');
                new_trigger_channel(1).digital_edge_polarity = fread(info_fid, 1, 'int16');
                new_channel(1).electrode_impedance_magnitude = fread(info_fid, 1, 'single');
                new_channel(1).electrode_impedance_phase = fread(info_fid, 1, 'single');
                
                if (channel_enabled)
                    switch (signal_type)
                        case 0
                            amplifier_channels(amplifier_index) = new_channel;
                            spike_triggers(amplifier_index) = new_trigger_channel;
                            amplifier_index = amplifier_index + 1;
                        case 1
                            aux_input_channels(aux_input_index) = new_channel;
                            aux_input_index = aux_input_index + 1;
                        case 2
                            supply_voltage_channels(supply_voltage_index) = new_channel;
                            supply_voltage_index = supply_voltage_index + 1;
                        case 3
                            board_adc_channels(board_adc_index) = new_channel;
                            board_adc_index = board_adc_index + 1;
                        case 4
                            board_dig_in_channels(board_dig_in_index) = new_channel;
                            board_dig_in_index = board_dig_in_index + 1;
                        case 5
                            board_dig_out_channels(board_dig_out_index) = new_channel;
                            board_dig_out_index = board_dig_out_index + 1;
                        otherwise
                            error('Unknown channel type');
                    end
                end
                
            end
        end
    end
    
    % Summarize contents of data file.
    num_amplifier_channels = amplifier_index - 1;
    num_aux_input_channels = aux_input_index - 1;
    num_supply_voltage_channels = supply_voltage_index - 1;
    num_board_adc_channels = board_adc_index - 1;
    num_board_dig_in_channels = board_dig_in_index - 1;
    num_board_dig_out_channels = board_dig_out_index - 1;
    
    fprintf(1, 'Found %d amplifier channel%s.\n', ...
        num_amplifier_channels, plural(num_amplifier_channels));
    fprintf(1, 'Found %d auxiliary input channel%s.\n', ...
        num_aux_input_channels, plural(num_aux_input_channels));
    fprintf(1, 'Found %d supply voltage channel%s.\n', ...
        num_supply_voltage_channels, plural(num_supply_voltage_channels));
    fprintf(1, 'Found %d board ADC channel%s.\n', ...
        num_board_adc_channels, plural(num_board_adc_channels));
    fprintf(1, 'Found %d board digital input channel%s.\n', ...
        num_board_dig_in_channels, plural(num_board_dig_in_channels));
    fprintf(1, 'Found %d board digital output channel%s.\n', ...
        num_board_dig_out_channels, plural(num_board_dig_out_channels));
    fprintf(1, 'Found %d temperature sensor channel%s.\n', ...
        num_temp_sensor_channels, plural(num_temp_sensor_channels));
    fprintf(1, '\n');
    
    % Determine how many samples the data file contains.
    
    % Each data block contains num_samples_per_data_block amplifier samples.
    bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_amplifier_channels;
    % Auxiliary inputs are sampled 4x slower than amplifiers
    bytes_per_block = bytes_per_block + (num_samples_per_data_block / 4) * 2 * num_aux_input_channels;
    % Supply voltage is sampled once per data block
    bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
    % Board analog inputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
    % Board digital inputs are sampled at same rate as amplifiers
    if (num_board_dig_in_channels > 0)
        bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
    end
    % Board digital outputs are sampled at same rate as amplifiers
    if (num_board_dig_out_channels > 0)
        bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
    end
    % Temp sensor is sampled once per data block
    if (num_temp_sensor_channels > 0)
       bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels; 
    end
 
    fclose(info_fid);
    
    %% Code written by Patrick to get "one file per channel" data into arrays

    %Figure out how much space to preallocate for the data, this is easily
    %done by checking the filesize of the time vector, and dividing by 4,
    %because there's 4 bytes in each time data point.
    time_id = fopen(folder_path + "\time.dat");
    fseek(time_id,0,"eof");
    file_size_of_time = ftell(time_id);
    num_samples_total = file_size_of_time/4;
    fseek(time_id,0,"bof");

    %Store the timestamp vector in an array
    time_stamp = fread(time_id,num_samples_total,"int32")';
    fclose(time_id);
    
    % Check for gaps in timestamps.
    num_gaps = sum(diff(time_stamp) ~= 1);
    if (num_gaps == 0)
        fprintf(1, 'No missing timestamps in data.\n');
    else
        fprintf(1, 'Warning: %d gaps in timestamp data found.  Time scale will not be uniform!\n', ...
            num_gaps);
    end

    %Store all the neuro data from the probe into an array
    if (isequal('ALL',desired_port))
        num_ports_desired = num_amplifier_channels;
    else
        num_ports_desired = 0;
        for check_ports = 1:num_amplifier_channels
            if (amplifier_channels(check_ports).port_prefix == desired_port)
                num_ports_desired = num_ports_desired + 1;
            end
        end
    end


    if (isequal('ALL',desired_port))
        desired_port_idx = 1:num_ports_desired;
    else
        desired_port_idx = find(vertcat(amplifier_channels.port_prefix) == desired_port);
    end

    amplifier_data = zeros(num_ports_desired,num_samples_total,'int16');
    %Each bit of data can be read on a different core for speed
    parfor iac = 1:length(desired_port_idx)        
        %In order for the data to be read, locate where it is, which is
        %done with folder_path, and knowing that each amplifier channel
        %begins with amp-, ends in .dat, and the middle names are stored in
        %amplifier_channels, which was generated earlier
        current_fid = fopen(folder_path+"\amp-"+amplifier_channels(desired_port_idx(iac)).native_channel_name+".dat");
        amplifier_data(iac,:) = fread(current_fid,num_samples_total,'int16') * 0.195;
        fclose(current_fid);
    end

    %Store all the auxillary input data
    %This is sampled at 1/4th the speed of amplifier channels
    %The vector is the same length, just there's multiple copies of the
    %same number 
    %For example, if amplifier is [1 2 3 1 2 3 2 4]
    %This one would be            [6 6 6 6 7 7 7 7]
    aux_input_data = zeros(num_aux_input_channels,num_samples_total,'uint16');
    for iaux = 1:num_aux_input_channels
        current_fid = fopen(folder_path + "\aux-" + aux_input_channels(iaux).native_channel_name+".dat");
        aux_input_data(iaux,:) = fread(current_fid,num_samples_total,'uint16') * 0.0000374;
        fclose(current_fid);
    end

    %Supply voltage
    %1/60th the speed of amplifier samples
    supply_voltage_data = zeros(num_supply_voltage_channels,num_samples_total,"uint16");
    parfor isv = 1:num_supply_voltage_channels
        current_fid = fopen(folder_path + "\vdd-" + supply_voltage_channels(isv).native_channel_name+".dat");
        supply_voltage_data(isv,:) = fread(current_fid,num_samples_total,'uint16') * 0.0000748;
        fclose(current_fid);
    end

    %Board ADC
    %Same rate as the amplifier data
    board_adc_data = zeros(num_board_adc_channels,num_samples_total,'single');
    parfor iadc = 1:num_board_adc_channels
        current_fid = fopen(folder_path + "\board-" + board_adc_channels(iadc).native_channel_name + ".dat");
        if (board_mode == 0)
            board_adc_data(iadc,:) = fread(current_fid,num_samples_total,'uint16') * 0.000050354;
        else
            board_adc_data(iadc,:) = (fread(current_fid,num_samples_total,'uint16')-32768) * 0.0003125;
        end
        fclose(current_fid);
    end

    %Board digital input data
    %Same rate as the amplifier
    board_dig_in_data = zeros(num_board_dig_in_channels,num_samples_total,'logical');
    for idin = 1:num_board_dig_in_channels
        current_fid = fopen(folder_path + "\board-" + board_dig_in_channels(idin).native_channel_name + ".dat");
        board_dig_in_data(idin,:) = fread(current_fid,num_samples_total,'uint16');
        fclose(current_fid);
    end

    %Board digital output data
    %Same rate as the amplifier
    board_dig_out_data = zeros(num_board_dig_out_channels,num_samples_total,'logical');
    parfor ido = 1:num_board_dig_out_channels
        current_fid = fopen(folder_path + "\board-" + board_dig_out_channels(ido).native_channel_name + ".dat");
        board_dig_out_data(ido,:) = fread(current_fid,num_samples_total,'uint16');
        fclose(current_fid);
    end

    %TODO FIGURE OUT SPIKE TRIGGERS

    %TODO There may or may not be temperature sensors, but if there are,
    %they should go here

    %Create a structure to output all the variables here
    probe_data.sampling_rate = sample_rate;
    probe_data.notes = notes;
    probe_data.frequency_parameters = frequency_parameters;

    if (data_file_main_version_number > 1)
        probe_data.reference_channel = reference_channel;
    end

    probe_data.amplifier_channels = amplifier_channels;
    probe_data.amplifier_data = amplifier_data;
    
    probe_data.time_stamp = time_stamp;
    probe_data.spike_triggers = spike_triggers;
    
    probe_data.aux_input_channels = aux_input_channels;
    probe_data.aux_input_data = aux_input_data;

    probe_data.supply_voltage_channels = supply_voltage_channels;
    probe_data.supply_voltage_data = supply_voltage_data;

    probe_data.board_adc_channels = board_adc_channels;
    probe_data.board_adc_data = board_adc_data;

    probe_data.board_dig_in_channels = board_dig_in_channels;
    probe_data.board_dig_in_data = board_dig_in_data;

    probe_data.board_dig_out_channels = board_dig_out_channels;
    probe_data.board_dig_out_data = board_dig_out_data;

    probe_data.selected_port = desired_port;
end