function task_data = VANDERBILT_PassiveGLOv2(task_data)

go = task_data.orientation(find(~isnan(task_data.orientation),1));
if      go == 135;  lo = 45;
elseif  go == 45;   lo = 135;
end

block_mat   = ["gloexp", "rndctl", "seqctl"];
seq_combos  = [ ...
    go, go, go, go; ...
    go, go, go, lo; ...
    go, go, lo, go; ...
    go, lo, go, go; ...
    lo, go, go, go; ...
    go, go, lo, lo; ...
    go, lo, go, lo; ...
    lo, go, go, lo; ...
    go, lo, lo, go; ...
    lo, go, lo, go; ...
    lo, lo, go, go; ...
    go, lo, lo, lo; ...
    lo, go, lo, lo; ...
    lo, lo, go, lo; ...
    lo, lo, lo, go; ...
    lo, lo, lo, lo ];

trial_start = find(contains(task_data.event_code_type, 'trial start'));
trial_end = [trial_start(2:end)-1; numel(task_data.start_time)];

temp_seq_type = {'gloexp'};
block_switch = 0;
for ii = 1 : numel(trial_start)

    presentations = find(contains(task_data.event_code_type(trial_start(ii):trial_end(ii)), 'task_event'));
    if numel(presentations) == 4
        complete_trial = 1;
        prez_seq = task_data.orientation(presentations+(trial_start(ii)-1))';
        task_data.seq_type(trial_start(ii):trial_end(ii)) = find(sum(seq_combos == prez_seq,2) == 4);
    else
        complete_trial = 0;
        task_data.seq_type(trial_start(ii):trial_end(ii)) = NaN;
    end

    task_data.presentation(trial_start(ii):trial_end(ii)) = NaN;
    if ~isempty(presentations)
        task_data.presentation(presentations+(trial_start(ii)-1)) = 1:numel(presentations);
    end

    if complete_trial
        if task_data.seq_type(trial_start(ii)) ~= 1 & ...
                task_data.seq_type(trial_start(ii)) ~= 2 & ...
                ~isnan(task_data.seq_type(trial_start(ii))) & ...
                ~block_switch
            block_switch = 1;
            temp_seq_type = {'rndctl'};
        end
    end

    task_data.sequence_type(trial_start(ii):trial_end(ii)) = ...
                repmat(temp_seq_type, numel(trial_start(ii):trial_end(ii)), 1);

    clear presentations prez_seq complete_trial
end

task_data.seq_type = task_data.seq_type';
task_data.sequence_type = task_data.sequence_type';
task_data.presentation = task_data.presentation';

for ii = numel(trial_start) : -1 : 1
    if task_data.seq_type(trial_start(ii)) ~= 1 & ...
            task_data.seq_type(trial_start(ii)) ~= 16 & ...
            ~isnan(task_data.seq_type(trial_start(ii)))

        task_data.sequence_type(trial_start(ii):end) = ...
            repmat({'seqctl'}, numel(trial_start(ii):numel(task_data.sequence_type)), 1);

        break

    end
end

task_data.gloexp = strcmp(task_data.sequence_type, 'gloexp');
task_data.rndctl = strcmp(task_data.sequence_type, 'rndctl');
task_data.seqctl = strcmp(task_data.sequence_type, 'seqctl');

% and code them relative to the combo matrix above
seq_go_type                        = 1;
seq_lo_type                        = 2;
seq_ilo_type                       = 15;
seq_igo_type                       = 16;

% identify useful sequences
task_data.go_seq                             = task_data.seq_type == seq_go_type;
task_data.lo_seq                             = task_data.seq_type == seq_lo_type;
task_data.igo_seq                            = task_data.seq_type == seq_igo_type;
task_data.ilo_seq                            = task_data.seq_type == seq_ilo_type;

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

if sum(strcmp(task_data.sequence_type, 'gloexp')) < sum(strcmp(task_data.sequence_type, 'rndctl'))
    warning('something seems to be wrong with sequence identification')
end

end