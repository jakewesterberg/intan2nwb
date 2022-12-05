function nwb = i2nDIO(pp, nwb, recdev)

% ADD FUNCTIONALITY TO GRAB DIGITAL CODES DYNAMICALLY. CURRENTLY RESTRICTED
% TO THE LAST 8 BITS OF THE DIGITAL OUTPUT.

digital_data = recdev.board_dig_in_data(end-7:end,:);
intan_code_times_unprocessed = find(sum(digital_data) > 0);
intan_code_times = nan(length(intan_code_times_unprocessed),1);
intan_code_values = nan(8,length(intan_code_times_unprocessed));

temp_ctr = 1;
intan_code_times(temp_ctr) = intan_code_times_unprocessed(temp_ctr);
intan_code_values(:,temp_ctr) = digital_data(:,intan_code_times(temp_ctr));
previous_value = intan_code_times(temp_ctr);
temp_ctr = temp_ctr + 1;
for jj = 2:length(intan_code_times_unprocessed)
    if ~(intan_code_times_unprocessed(jj) == previous_value + 1)
        intan_code_times(temp_ctr) = intan_code_times_unprocessed(jj);
        intan_code_values(:,temp_ctr) = digital_data(:,intan_code_times(temp_ctr)+1);
        temp_ctr = temp_ctr + 1;
    end
    previous_value = intan_code_times_unprocessed(jj);
end

intan_code_times = intan_code_times(1:temp_ctr-1) ./ recdev.sampling_rate;
intan_code_values = intan_code_values(:,1:temp_ctr-1);
intan_code_values = bit2int(flip(intan_code_values),8)';

% Dientangle event codes...
event_data = identEvents(intan_code_values, intan_code_times);

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
    nwb.intervals.set(event_data{jj}.task, trials); clear trials
end

%nwbExport(nwb, [pp.NWB_DATA nwb.identifier '.nwb']);

end