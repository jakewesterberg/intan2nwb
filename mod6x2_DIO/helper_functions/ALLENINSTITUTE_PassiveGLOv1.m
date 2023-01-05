function task_data = ALLENINSTITUTE_PassiveGLOv1(nwb)

task_data.task = 'passive_glo';

% GLO presentation times
task_data.start_time                                 = nwb.intervals.get('init_grating_presentations').start_time.data(:);
task_data.stop_time                                  = nwb.intervals.get('init_grating_presentations').stop_time.data(:);

% GLO stimulus information
task_data.orientation                                = nwb.intervals.get('init_grating_presentations').vectordata.get('orientation').data(:); % orientation
task_data.drift_rate                                 = nwb.intervals.get('init_grating_presentations').vectordata.get('temporal_frequency').data(:); % temporal frequency (drift rate)  
task_data.spatial_freq                               = nwb.intervals.get('init_grating_presentations').vectordata.get('spatial_frequency').data(:); % spatial frequency
task_data.contrast                                   = nwb.intervals.get('init_grating_presentations').vectordata.get('contrast').data(:); % stimulus contrast

% first presentation is always a local oddball trial
go              = task_data.orientation(1);
lo              = task_data.orientation(4);

% all possible sequence combinations
seq_combos                         = [ go, go, go, go; go, go, go, lo; ...
                                       go, go, lo, go; go, lo, go, go; ...
                                       lo, go, go, go; go, go, lo, lo; ...
                                       go, lo, go, lo; lo, go, go, lo; ...
                                       go, lo, lo, go; lo, go, lo, go; ...
                                       lo, lo, go, go; go, lo, lo, lo; ...
                                       lo, go, lo, lo; lo, lo, go, lo; ...
                                       lo, lo, lo, go; lo, lo, lo, lo ];

% count total trials in glo
total_trials                       = numel(task_data.start_time);

% identify important sequence types
seq_go                             = [go, go, go, go];
seq_igo                            = [lo, lo, lo, lo];
seq_lo                             = [go, go, go, lo];

% and code them relative to the combo matrix above
seq_go_type                        = 1;
seq_lo_type                        = 2;
seq_ilo_type                       = 15;
seq_igo_type                       = 16;

% find glo sequence info
task_data.seq_type                           = [];
for i = 1 : 5 : total_trials
    if (~floor(sum(task_data.orientation(i:i+3)==seq_go')/4) & ...
            ~floor(sum(task_data.orientation(i:i+3)==seq_lo')/4)) & ...
            ~isfield(task_data, 'gloexp')

        % index the positions of main experiment vs. control sequence presentations
        task_data.gloexp                     = zeros(total_trials,1); task_data.gloexp(1:i-1) = 1; % main glo sequences
        rndctl_start = i;

    end

    % code the sequence types in case we determine other seq combos are
    % interesting
    task_data.seq_type                       = [task_data.seq_type; repmat(find(sum(seq_combos == task_data.orientation(i:i+3)',2)==4),5,1)]; % give each sequence a code relative combo matrix

end

detected_start_of_late_glo = 0;
for i = total_trials-4 : -5 : 1
    if (~floor(sum(task_data.orientation(i:i+3)==seq_go')/4) & ...
            ~floor(sum(task_data.orientation(i:i+3)==seq_lo')/4)) & ...
            ~detected_start_of_late_glo

        % index the positions of main experiment vs. control sequence presentations
        if i ~= total_trials-4 & i ~= total_trials-9
            task_data.gloexp(i+5:end)            = 1; % main glo sequences
            seqctl_end = i+4;
        else
            seqctl_end = total_trials;
        end
        detected_start_of_late_glo = 1;

    end

    if (~floor(sum(task_data.orientation(i:i+3)==seq_go')/4) & ...
            ~floor(sum(task_data.orientation(i:i+3)==seq_igo')/4)) & ...
            detected_start_of_late_glo

          task_data.rndctl                     = zeros(total_trials,1); task_data.rndctl(rndctl_start:i+4) = 1; % randomized control sequences
          task_data.seqctl                     = zeros(total_trials,1); task_data.seqctl(i+5:seqctl_end) = 1; % sequenced control with alternating [g g g g, l l l l, etc.]

        break

    end
end

% identify useful sequences
task_data.go_seq                             = task_data.seq_type == seq_go_type;
task_data.lo_seq                             = task_data.seq_type == seq_lo_type;
task_data.igo_seq                            = task_data.seq_type == seq_igo_type;
task_data.ilo_seq                            = task_data.seq_type == seq_ilo_type;

% index presentation numbers within sequence
task_data.presentation                       = zeros(total_trials,1); 
task_data.presentation(1:5:end)              = 1; % first presentation in a sequence (no adaptation)
task_data.presentation(2:5:end)              = 2;               
task_data.presentation(3:5:end)              = 3; % stimulus prior to the oddball              
task_data.presentation(4:5:end)              = 4; % sequence position we are most interested in (where the oddball occurs)          

task_data.trial_num                          = sort(repmat((1:total_trials/5)', 5, 1));

% predetermine some useful combinations
task_data.go_gloexp                          = task_data.seq_type == seq_go_type & task_data.gloexp & task_data.presentation==4; % global oddball presentations
task_data.lo_gloexp                          = task_data.seq_type == seq_lo_type & task_data.gloexp & task_data.presentation==4; % local oddball presentations

task_data.go_rndctl                          = task_data.seq_type == seq_go_type & task_data.rndctl & task_data.presentation==4; % 'global oddball' presentation in random control
task_data.lo_rndctl                          = task_data.seq_type == seq_lo_type & task_data.rndctl & task_data.presentation==4; % 'local oddball' presentation in random control
task_data.igo_rndctl                         = task_data.seq_type == seq_igo_type & task_data.rndctl & task_data.presentation==4; % inverse 'global oddball' presentation in random control [l l l g] insead of [g g g l]
task_data.ilo_rndctl                         = task_data.seq_type == seq_ilo_type & task_data.rndctl & task_data.presentation==4; % inverse 'local oddball' presentation in random control [l l l l] insead of [g g g g]

task_data.go_seqctl                          = task_data.seq_type == seq_go_type & task_data.seqctl & task_data.presentation==4; % 'global oddball' presentation in sequence control
task_data.igo_seqctl                         = task_data.seq_type == seq_igo_type & task_data.seqctl & task_data.presentation==4; % inverse 'local oddball' presentation in sequence control [l l l l] insead of [g g g g]

% Convert to logicals
task_data.go_seq                             = logical(task_data.go_seq);
task_data.lo_seq                             = logical(task_data.lo_seq);
task_data.igo_seq                            = logical(task_data.igo_seq);
task_data.ilo_seq                            = logical(task_data.ilo_seq);
task_data.gloexp                             = logical(task_data.gloexp);
task_data.rndctl                             = logical(task_data.rndctl);
task_data.seqctl                             = logical(task_data.seqctl);

task_data.early_glo                          = task_data.gloexp; task_data.early_glo(find(task_data.rndctl,1,"first"):end) = 0;
task_data.late_glo                           = task_data.gloexp; task_data.late_glo(1:find(task_data.rndctl,1,"first")) = 0;

end