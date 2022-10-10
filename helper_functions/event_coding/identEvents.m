function event_data = identEvents(codes, times)
event_data = {};

task_ident_code     = 100;
unique_task_idents  = { ...
    ...
    'passive_glo',          1,      ...
    'active_glo',           2,      ...
    'attention_glo',        3,      ...
    ...
    };

tasks_completed = unique(codes(find(codes==100)+1))';

if isempty(tasks_completed)

    event_data{1} = task_passive_glo_v1(codes, times);
    
else

    for ii = tasks_completed
        switch ii
            case 1 % passive_glo
            case 2 % active_glo
            case 3 % attention_glo
        end
    end
end
end