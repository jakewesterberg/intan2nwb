function task_data = VANDERBILT_BastosGeneralizedEventDecode(trialified, task_input)

ID = '1Mn6JhVGvZK30M_1pl5i4BrDiwX4rGLaJ0lYybH7EFn0';
url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s', ID);
event_info = webread(url_name);

ii_ctr = 1;
for ii = 1 : numel(trialified)

    primary_codes = trialified{ii}.codes(1:find(trialified{ii}.codes == 84, 1, "first")-1);
    primary_times = trialified{ii}.times(1:find(trialified{ii}.codes == 84, 1, "first")-1);

    secondary_codes = trialified{ii}.codes(find(trialified{ii}.codes == 1001, 1, "first"):end-1);
    secondary_codes = [secondary_codes(1:2:end-1), secondary_codes(2:2:end)];
    stimulus_specific_codes = secondary_codes(find(secondary_codes(:,1)==2001):end,:);
    secondary_codes = secondary_codes(1:find(secondary_codes(:,1)==2001)-1,:);

    stimuli = find(stimulus_specific_codes(:,1) == 2001);
    n_stimuli = numel(stimuli);
    if n_stimuli > 0
        for jj = 1:n_stimuli
            if jj < n_stimuli
                ssct{jj} = stimulus_specific_codes(stimuli(jj):stimuli(jj+1)-1,:);
            else
                ssct{jj} = stimulus_specific_codes(stimuli(jj):end,:);
            end
        end
        stimulus_specific_codes = ssct; clear ssct;
    end

    for jj = 1 : size(secondary_codes, 1)

        if secondary_codes(jj,1) > 1000 % bug in output code
            sci = find(event_info.SECONDARY_EVENT_CODE == secondary_codes(jj,1));
        else
            if abs(secondary_codes(jj-1,1) - secondary_codes(jj+1,1)) == 2
                sci = find(event_info.SECONDARY_EVENT_CODE == secondary_codes(jj-1,1)) + 1;
            else
                warning('cannot parse secondary event code...pairings misaligned, trying my best')
                continue
            end
        end
        scs.(lower(event_info.SECONDARY_EVENT_CODE_DESCRIPTION{sci})) = secondary_codes(jj,2);
    end
    scs_fields = fields(scs);

    for jj = 1 : numel(scs_fields)
        if strcmp(scs_fields{jj}, 'task')
            tci = find(event_info.TASK_SPECIFIC_EVENT_CODE == scs.(scs_fields{jj}));
            if isempty(tci)
                tci = find(strcmp(lower(event_info.TASK_SPECIFIC_EVENT_CODE_DESCRIPTION), scs.(scs_fields{jj})));
            end
            if isempty(tci)
                scs.(scs_fields{jj}) = 'TASK_CODE_CORRUPTED';
            else
                scs.(scs_fields{jj}) = lower(event_info.TASK_SPECIFIC_EVENT_CODE_DESCRIPTION{tci});
                break
            end
        end
    end

    time_on_code = event_info.STIMULUS_SPECIFIC_EVENT_CODE( ...
        find(strcmp(lower(event_info.STIMULUS_SPECIFIC_EVENT_CODE_DESCRIPTION), "start_time")));

    for jj = 1 : numel(primary_codes)

        if (primary_codes(jj) < 100 | primary_codes(jj) > 250) & primary_codes(jj) ~= 10

            if primary_codes(jj) == 18
                task_data.stop_time(track_start_idx) = primary_times(jj);
                clear track_start_idx
            else

                task_data.codes(ii_ctr) = primary_codes(jj);
                task_data.start_time(ii_ctr) = primary_times(jj);
                task_data.stop_time(ii_ctr) = primary_times(jj);
                task_data.trial_num(ii_ctr) = ii;

                task_data.event_code_type{ii_ctr} = uncell(event_info.PRIMARY_EVENT_CODE_DESCRIPTION( ...
                    find(event_info.PRIMARY_EVENT_CODE == primary_codes(jj))),1);

                for kk = 1 : numel(scs_fields)
                    if isnumeric(scs.(scs_fields{kk}))
                        task_data.(scs_fields{kk})(ii_ctr) = scs.(scs_fields{kk});
                    else
                        task_data.(scs_fields{kk}){ii_ctr} = scs.(scs_fields{kk});
                    end
                end

                if primary_codes(jj) == 9
                    track_start_idx = ii_ctr;
                end

            end

            ii_ctr = ii_ctr + 1;

        end
    end

    for kk = 1 : numel(stimulus_specific_codes)

        if size(stimulus_specific_codes{kk}, 1) == 1
            continue
        end

        event_appeared = 1;
        sscs = [];

        for mm = 1 : size(stimulus_specific_codes{kk}, 1)
            if stimulus_specific_codes{kk}(mm,1) > 2000 % bug in output code
                ssci = find(event_info.STIMULUS_SPECIFIC_EVENT_CODE == stimulus_specific_codes{kk}(mm,1));
            else
                if abs(stimulus_specific_codes{kk}(mm-1,1) - stimulus_specific_codes{kk}(mm+1,1)) == 2
                    ssci = find(event_info.STIMULUS_SPECIFIC_EVENT_CODE == stimulus_specific_codes{kk}(mm-1,1)) + 1;
                else
                    warning('cannot parse stimulus specific event code...pairings misaligned, trying my best')
                    continue
                end
            end
            sscs.(lower(event_info.STIMULUS_SPECIFIC_EVENT_CODE_DESCRIPTION{ssci})) = ...
                stimulus_specific_codes{kk}(mm,2) * event_info.STIMULUS_SPECIFIC_EVENT_CODE_MULTIPLIER(ssci);
        end
        sscs_fields = fields(sscs);

        for mm = 1 : numel(sscs_fields)
            if event_appeared
                if isnumeric(sscs.(sscs_fields{mm}))
                    task_data.(sscs_fields{mm})(ii_ctr) = sscs.(sscs_fields{mm});
                    data_available_tracker.(sscs_fields{mm})(ii_ctr) = 1;
                else
                    task_data.(sscs_fields{mm}){ii_ctr} = sscs.(sscs_fields{mm});
                    data_available_tracker.(sscs_fields{mm})(ii_ctr) = 1;
                end

                if strcmp(sscs_fields{mm}, 'start_time')
                    try
                        task_data.start_time(ii_ctr) = trialified{ii}.times( ...
                            find(trialified{ii}.codes == task_data.start_time(ii_ctr), 1, "first"));
                    catch
                        % event never came to pass
                        event_appeared = 0;
                    end
                end

                if strcmp(sscs_fields{mm}, 'stop_time')
                    task_data.stop_time(ii_ctr) = trialified{ii}.times( ...
                        find(trialified{ii}.codes == task_data.stop_time(ii_ctr), 1, "first"));
                end
            end
        end

        if event_appeared
            if numel(task_data.start_time) < ii_ctr
                %warning('missing a start time! on early sessions this is okay due to the fixation code bug...')

                find_fix = find(trialified{ii}.codes == 10, 1, "first");
                if isempty(find_fix); find_fix = find(trialified{ii}.codes == 100, 1, "first"); end
                if isempty(find_fix); find_fix = find(trialified{ii}.codes == 101, 1, "first"); end

                task_data.start_time(ii_ctr) = trialified{ii}.times(find_fix);
                task_data.stop_time(ii_ctr) = trialified{ii}.times(find_fix);

            end

            task_data.codes(ii_ctr) = 99 + kk;
            task_data.trial_num(ii_ctr) = ii;

            task_data.event_code_type{ii_ctr} = lower(uncell(event_info.PRIMARY_EVENT_CODE_DESCRIPTION( ...
                find(event_info.PRIMARY_EVENT_CODE == task_data.codes(ii_ctr))),1));

            for mm = 1 : numel(scs_fields)
                if isnumeric(scs.(scs_fields{mm}))
                    task_data.(scs_fields{mm})(ii_ctr) = scs.(scs_fields{mm});
                else
                    task_data.(scs_fields{mm}){ii_ctr} = scs.(scs_fields{mm});
                end
            end

            if isfield(task_data, 'is_fixation_dot')
                if numel(task_data.is_fixation_dot) == ii_ctr
                    if task_data.is_fixation_dot(ii_ctr) == 1
                        task_data.event_code_type{ii_ctr} = 'fix cue appearance';
                    end
                end
            end

            ii_ctr = ii_ctr + 1;

        end
    end
end

task_data_fields = fields(data_available_tracker);
for ii = 1 : numel(task_data_fields)
    if ~strcmp(task_data_fields{ii}, 'start_time') & ~strcmp(task_data_fields{ii}, 'stop_time')
        task_data.(task_data_fields{ii})(~data_available_tracker.(task_data_fields{ii})) = NaN;
    end
    if numel(task_data.(task_data_fields{ii})) ~= numel(task_data.codes)
        task_data.(task_data_fields{ii}) = [task_data.(task_data_fields{ii}) ...
            nan(1, numel(task_data.codes) - numel(task_data.(task_data_fields{ii})))];
    end
end

task_data_fields = fields(task_data);
[~, sort_idx] = sort(task_data.start_time);
for ii = 1 : numel(task_data_fields)
    if size(task_data.(task_data_fields{ii}), 2) > 1
        task_data.(task_data_fields{ii}) = task_data.(task_data_fields{ii})';
    end
%     task_data.(task_data_fields{ii}) = task_data.(task_data_fields{ii})(sort_idx);
end

if exist('task_input', 'var')
    if isnumeric(task_input)
        task_input = lower(event_info.TASK_SPECIFIC_EVENT_CODE_DESCRIPTION{task_input});
    end
    switch task_input
        case 'passive_glo'
            task_data = VANDERBILT_PassiveGLOv2(task_data);
    end
end

end