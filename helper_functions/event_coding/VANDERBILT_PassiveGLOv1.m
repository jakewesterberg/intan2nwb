function task_data = taskPassiveGLOv1(codes, times)

ii_ctr = 1;

block_mat = ["gloexp", "rndctl", "seqctl"];
condition_mat = [ ...
    1 45 45 45 135; ...
    1 45 45 45 135; ...
    1 45 45 45 45; ...
    2 45 45 45 45; ...
    2 45 45 45 45; ...
    2 45 45 45 45; ...
    2 135 135 135 135; ...
    2 135 135 135 135; ...
    2 135 135 135 135; ...
    2 45 45 45 135; ...
    2 45 45 45 135; ...
    2 45 45 45 135; ...
    2 135 135 135 45; ...
    2 135 135 135 45; ...
    2 135 135 135 45; ...
    2 45 45 135 135; ...
    2 45 135 135 135; ...
    2 135 135 45 45; ...
    2 135 45 45 45; ...
    2 135 45 45 135; ...
    2 45 135 135 45; ...
    2 45 135 45 45; ...
    2 45 45 135 45; ...
    2 45 135 45 135; ...
    2 135 45 135 45; ...
    2 135 45 135 135; ...
    2 135 135 45 135; ...
    3 135 135 135 135; ...
    3 45 45 45 45];

trial_ct = 0;
for ii = 1 : numel(codes)

    if codes(ii)>100 & codes(ii)<150
        current_trial_seq = condition_mat(codes(ii)-100, 1);
        current_trial_type = condition_mat(codes(ii)-100, 2:5);
        temp_note = codes(ii) - 100;
        trial_ct = trial_ct + 1;
        continue
    end

    if codes(ii) == 10 | codes(ii) == 255 | codes(ii) == 20 | codes(ii) == 22 | codes(ii) == 24 | codes(ii) == 30 | codes(ii) == 40
        task_data.codes(ii_ctr)         = codes(ii);
        task_data.start_time(ii_ctr)   = times(ii);
        task_data.trial_num(ii_ctr)     = trial_ct;
        try
            task_data.stop_time(ii_ctr)    = times(ii+1);
        catch
            task_data.stop_time(ii_ctr)    = times(ii);
        end
        switch codes(ii)
            case 10
                task_data.orientation(ii_ctr) = NaN;
                task_data.presentation(ii_ctr) = NaN;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "fix cue appearance";
            case 255
                task_data.orientation(ii_ctr) = NaN;
                task_data.presentation(ii_ctr) = NaN;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "fixation made";
            case 20
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(1));
                task_data.presentation(ii_ctr) = 1;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "presentation 1";
            case 22
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(2));
                task_data.presentation(ii_ctr) = 2;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "presentation 2";
            case 24
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(3));
                task_data.presentation(ii_ctr) = 3;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "presentation 3";
            case 30
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(4));
                task_data.presentation(ii_ctr) = 4;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "presentation 4";
            case 40
                task_data.orientation(ii_ctr) = NaN;
                task_data.presentation(ii_ctr) = NaN;
                task_data.sequence_type{ii_ctr} = block_mat(current_trial_seq);
                task_data.notes{ii_ctr} = temp_note;
                task_data.event_code_type{ii_ctr} = "reward";
        end
        ii_ctr = ii_ctr +1;
    end
end

task_data.correct = zeros(numel(task_data.codes), 1, 'int8');
u_trials = unique(task_data.trial_num);
for ii = u_trials
    temp_inds = sum(task_data.codes == 40 & task_data.trial_num == ii);
    if temp_inds
        task_data.correct(task_data.trial_num == ii) = int8(1);
    end
end

task_data.codes             = task_data.codes';
task_data.trial_num         = task_data.trial_num';
task_data.presentation      = task_data.presentation';
task_data.start_time        = task_data.start_time';
task_data.stop_time         = task_data.stop_time';
task_data.orientation       = task_data.orientation';
task_data.sequence_type     = task_data.sequence_type';
task_data.notes             = task_data.notes';
task_data.event_code_type   = task_data.event_code_type';

task_data.task = 'passive_glo';

end