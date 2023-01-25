function event_data = identEvents(codes, times, task_ident_code, task_trial_start_code, task_trial_end_code)

%Compiling function that directs the detected event codes specifiers to the
%correct event coding function. event coding function should output a
%structure with  nx1 fields with the names that will be assigned to the nwb
%field with the .task name.

event_data = {};

if nargin < 3
    task_ident_code         = 1001;
    task_trial_start_code   = 9;
    task_trial_end_code     = [85; 85; 85];
end

% missing trial start code bug fix
for ii = 2 : numel(codes)-1
    if codes(ii) == 18 & codes(ii-1) == 18 & codes(ii+1) ~= 9 & codes(ii+1) ~= 18
        codes(ii) = 9;
    end
end

codes_plus_2 = [codes(3:end); NaN; NaN];

tasks_completed = unique(codes(find(codes==task_ident_code & codes_plus_2==task_ident_code+1)+1))';

end_codes = zeros(1, numel(codes),'logical');
for i = 1 : numel(codes)-numel(task_trial_end_code)+1
    if sum(codes(i:i+numel(task_trial_end_code)-1) == task_trial_end_code) == numel(task_trial_end_code)
        end_codes(i) = 1;
    end
end

trialified = [];
trials = find(codes==task_ident_code & codes_plus_2==task_ident_code+1);
assigned_task = zeros(1,numel(trials));
for i = 1 : numel(trials)
    if i == 1
        trial_begin = trials(i)-1;
    else
        trial_begin = min(abs(trials(i) - find(codes(1:trials(i))==task_trial_start_code)));
    end
    trial_end = find(end_codes(trials(i):end), 1, "first")-2;

    assigned_task(i) = codes(trials(i)+1);
    if i < numel(trials) & i > 1
        if assigned_task(i) ~= assigned_task(i-1) & ...
                assigned_task(i) ~= codes(trials(i+1)+1) & ...
                assigned_task(i-1) == codes(trials(i+1)+1)
            assigned_task(i) = assigned_task(i-1);
        end
    end

    try
        trialified{i}.codes = codes(trials(i) - trial_begin:trials(i) + trial_end);
        trialified{i}.times = times(trials(i) - trial_begin:trials(i) + trial_end);
    catch
        trialified{i}.codes = codes(trials(i) - trial_begin:end);
        trialified{i}.times = times(trials(i) - trial_begin:end);
    end

    % correct the early-days-early-event bug
    if trialified{i}.codes(2) == 102
        trialified{i}.codes = [trialified{i}.codes(1); trialified{i}.codes(3:end)];
        trialified{i}.times = [trialified{i}.times(1); trialified{i}.times(3:end)];
    end
    if trialified{i}.codes(1) == 9 & trialified{i}.codes(2) >= 100 & trialified{i}.codes(3) == 10
        trialified{i}.codes = [trialified{i}.codes(1); trialified{i}.codes(3:end)];
        trialified{i}.times = [trialified{i}.times(1); trialified{i}.times(3:end)];
    end

end

% fix strobe position bug
if trialified{1}.codes(1) == 18 & trialified{1}.codes(2) == 9
    trialified{1}.codes = trialified{1}.codes(2:end);
    trialified{1}.times = trialified{1}.times(2:end);
end
if trialified{1}.codes(1) == 18 & trialified{1}.codes(2) == 18
    trialified{1}.codes(2) = 9;
    trialified{1}.codes = trialified{1}.codes(2:end);
    trialified{1}.times = trialified{1}.times(2:end);
end
if trialified{1}.codes(1) == 18 & (trialified{1}.codes(2) ~= 18 & trialified{1}.codes(2) ~= 9)
    trialified{1}.codes(1) = 9;
end
if trialified{1}.codes(2) == 101 & trialified{1}.codes(3) == 10
    trialified{1}.codes = [trialified{1}.codes(1); trialified{1}.codes(3:end)];
    trialified{1}.times = [trialified{1}.times(1); trialified{1}.times(3:end)];
end

tasks = unique(assigned_task);
for ii = 1:numel(tasks)
    event_data{ii} = VANDERBILT_BastosGeneralizedEventDecode( ...
        trialified(assigned_task == tasks(ii)), ...
        tasks(ii));
end
end