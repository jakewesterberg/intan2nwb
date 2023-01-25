function nwb2 = i2nAIC_events_only(pp, nwb, recording_info, ii)

nwb2                                 = NwbFile;

% reformat existing
nwb2.identifier                      = recording_info.Identifier{ii};
nwb2.general_experimenter            = recording_info.Investigator{ii};
nwb2.general_session_id              = recording_info.Identifier{ii};
nwb2.general_experiment_description  = recording_info.Experiment_Description{ii};

nwb2.general_extracellular_ephys_electrodes = nwb.general_extracellular_ephys_electrodes;

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