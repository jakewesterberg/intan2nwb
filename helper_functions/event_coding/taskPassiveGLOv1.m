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

for ii = 1 : numel(codes)

    if codes(ii)>100
        current_trial_seq = condition_mat(codes(ii)-100, 1);
        current_trial_type = condition_mat(codes(ii)-100, 2:5);
        temp_note = codes(ii) - 100;
        continue
    end

    task_data.codes(ii_ctr)         = codes(ii);
    task_data.start_times(ii_ctr)   = times(ii);
    if codes(ii) == 10 | codes(ii) == 255 | codes(ii) == 20 | codes(ii) == 22 | codes(ii) == 24 | codes(ii) == 30 | codes(ii) == 40
        try
            task_data.stop_times(ii_ctr)    = times(ii+1);
        catch
            task_data.stop_times(ii_ctr)    = times(ii);
        end
        switch codes(ii)
            case 10
                task_data.orientation(ii_ctr) = NaN;
                task_data.sequence_type = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'fix cue appearance';
            case 255
                task_data.orientation(ii_ctr) = NaN;
                task_data.sequence_type = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'fixation made';
            case 20
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(1));
                task_data.sequence_type(ii_ctr) = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'presentation 1';
            case 22
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(2));
                task_data.sequence_type(ii_ctr) = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'presentation 2';
            case 24
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(3));
                task_data.sequence_type(ii_ctr) = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'presentation 3';
            case 30
                task_data.orientation(ii_ctr) = condition_mat(current_trial_type(4));
                task_data.sequence_type(ii_ctr) = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'presentation 4';
            case 40
                task_data.orientation(ii_ctr) = NaN;
                task_data.sequence_type(ii_ctr) = block_mat(current_trial_seq);
                task_data.notes(ii_ctr) = temp_note;
                task_data.description(ii_ctr) = 'reward';
        end
    end
    ii_ctr = ii_ctr +1;
end

task_data.task = 'passive_glo';

end