function opto_info = ALLENINSTITUTE_Optotaggingv1(nwb)

% optotagging events
opto_info.start_time        = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').start_time.data(:);
opto_info.stop_time         = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').stop_time.data(:);

% stimulation info
opto_info.level             = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('level').data(:);
opto_info.type              = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('stimulus_name').data(:);

% count opto total trials
opto_info.presentation      = (1:numel(opto_info.start_time))';

opto_info.task              = 'optotagging';

end