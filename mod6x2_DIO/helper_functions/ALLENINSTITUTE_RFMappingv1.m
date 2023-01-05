function rf = ALLENINSTITUTE_RFMappingv1(nwb)

% Jake Westerberg
% Vanderbilt University
% jakewesterberg@gmail.com

rf.task =  'rf_mapping';

% grab rf trial times
rf.start_time               = nwb.intervals.get('create_receptive_field_mapping_presentations').start_time.data(:);
rf.stop_time                = nwb.intervals.get('create_receptive_field_mapping_presentations').stop_time.data(:);

% grab rf trial stim info
rf.x_pos                    = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('x_position').data(:);
rf.y_pos                    = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('y_position').data(:);
rf.size                     = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('size').data(:);
rf.phase                    = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('phase').data(:);
rf.orientation              = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('orientation').data(:);
rf.drift_rate               = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('temporal_frequency').data(:);
rf.spatial_freq             = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('spatial_frequency').data(:);
rf.contrast                 = nwb.intervals.get('create_receptive_field_mapping_presentations').vectordata.get('contrast').data(:);

% count total rf trials
rf.presentation             = (1:numel(rf.start_time))';

end