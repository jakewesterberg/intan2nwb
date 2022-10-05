%How many milliseconds before to have, just in case
BUFFER = 500;

%[ml_data,config] = mlconcatenate('220909_Leo_glo_controlA.bhv2');
%load('LEO_glo_v_probe_data_0909_2022.mat')
% folder_intan = dir('E:\~Jake\Glo VProbe\GLO_CONTROL_VPROBE_220910_220909_104721');
% intan_data = read_intan_data(folder_intan(1).folder,'Z');
% 
% digital_data = intan_data.board_dig_in_data;
% first_digital_signal = find(sum(digital_data) > 0);
% first_digital_signal = first_digital_signal(1);
% 
% total_codes = 0;
% ml_code_times_concat = [];
% ml_code_values_concat = [];
% for ii = 1:length(ml_data.BehavioralCodes)
%     total_codes = total_codes + length(ml_data.BehavioralCodes(ii).CodeTimes);
%     ml_code_times_concat = cat(1,ml_code_times_concat,ml_data.BehavioralCodes(ii).CodeTimes);
%     ml_code_values_concat = cat(1,ml_code_values_concat,ml_data.BehavioralCodes(ii).CodeNumbers);
% end
% 
% digital_data = intan_data.board_dig_in_data;
% intan_code_times_unprocessed = find(sum(digital_data) > 0);
% intan_code_times = nan(length(intan_code_times_unprocessed),1);
% intan_code_values = nan(8,length(intan_code_times_unprocessed));
% 
% counter = 1;
% intan_code_times(counter) = intan_code_times_unprocessed(counter);
% intan_code_values(:,counter) = digital_data(:,intan_code_times(counter));
% previous_value = intan_code_times(counter);
% counter = counter + 1;
% for ii = 2:length(intan_code_times_unprocessed)
%     if(intan_code_times_unprocessed(ii) == previous_value + 1)
%         %Do nothing
%     else
%         intan_code_times(counter) = intan_code_times_unprocessed(ii);
%         intan_code_values(:,counter) = digital_data(:,intan_code_times(counter)+1);
%         counter = counter + 1;
%     end
%     previous_value = intan_code_times_unprocessed(ii);
% end

intan_code_times = intan_code_times(1:counter-1);
intan_code_values = intan_code_values(:,1:counter-1);

%Convert from binary into decimal
intan_code_values = bit2int(flip(intan_code_values),8)';

total_codes = 0;
ml_code_times_concat = [];
ml_code_values_concat = [];
for ii = 1:length(ml_data.BehavioralCodes)
    total_codes = total_codes + length(ml_data.BehavioralCodes(ii).CodeTimes);
    ml_code_times_concat = cat(1,ml_code_times_concat,ml_data.BehavioralCodes(ii).CodeTimes);
    ml_code_values_concat = cat(1,ml_code_values_concat,ml_data.BehavioralCodes(ii).CodeNumbers);
end

%Check for problems
if (total_codes ~= length(ml_code_values_concat))
    fprintf('WARNING: DATA INCONSISTENT, (EVENTCODES DIFFERENT)')
    return;
end

if (sum(ml_code_values_concat - intan_code_values)~=0)
    fprintf('WARNING: DATA INCONSISTENT, (EVENTCODES DIFFERENT)');
    return
end

timing_difference = mean(abs(diff(ml_code_times_concat) - diff(intan_code_times/((intan_data.sampling_rate/1000)))));

if timing_difference > 0.1
    sprintf('WARNING: TIMINGS INCONSISTENT %d\n',timing_difference);
end

%Begin Processing
%Downsample pupil data
eye_data_all = downsample(intan_data.board_adc_data',30)';
new_timestamps = round(intan_code_times/30);

%Figure out each trial
num_trials = sum(intan_code_values == 9);

%EventCodes
event_codes = ml_data.BehavioralCodes;

%Create place to put the data
event_code_time_stamps = nan(num_trials,20);
eye_data_per_trial = nan(3,20000,num_trials);
low_pass_data_per_trial = nan(20000,32, num_trials);
spike_data_per_trial = nan(20000,32, num_trials);
eye_data_time_stamps = nan(num_trials,20);

%Fill timestamps
for ii = 1:num_trials
    current_times = round(event_codes(ii).CodeTimes);
    current_times = current_times - current_times(1) + BUFFER;
    event_code_time_stamps(ii,1:length(current_times)) = current_times;
end



%Start filling the pupil data
end_loc = 1;
for ii = 1:num_trials
    found_end = false;
    start_timestamp = new_timestamps(end_loc) - BUFFER;
    while(~found_end)
        end_loc = end_loc + 1;
        %Signifies the end of a trial
        if (intan_code_values(end_loc) == 18)
            end_timestamp = new_timestamps(end_loc);
            end_loc = end_loc + 1;
            found_end = true;
        end
    end
    %fill in the eye data for this trial
    difference_in_timing = end_timestamp - start_timestamp;
    eye_data_per_trial(:,1:difference_in_timing+1,ii) = ...
        eye_data_all(:,start_timestamp:end_timestamp);
    low_pass_data_per_trial(1:difference_in_timing+1,:,ii) = ...
        low_pass_data(:,start_timestamp:end_timestamp)';
    spike_data_per_trial(1:difference_in_timing+1,:,ii) = ...
        spike_data(:,start_timestamp:end_timestamp)';
end

%Fill timestamps
for ii = 1:num_trials
    current_times = round(event_codes(ii).CodeTimes);
    current_times = current_times - current_times(1) + BUFFER;
    eye_data_time_stamps(ii,1:length(current_times)) = current_times;
end



correct_trials = ml_data.TrialError == 0;

%Local Or Globall oddball
condition = ml_data.Condition;

%Create file title
subject_name = config.SubjectName;

date_part = char(strjoin(string(ml_data.TrialDateTime(1,1:5)),'\\'));

date_part = strrep(date_part,'\','_');

%save(string(subject_name)+"JAKE_EPOCH_VPROBE"+string(date_part),"event_code_time_stamps","eye_data_per_trial","low_pass_data_per_trial","spike_data_per_trial","event_codes","correct_trials","condition","-v7.3")
