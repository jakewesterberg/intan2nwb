function event_data = identEvents(codes, times)

%Compiling function that directs the detected event codes specifiers to the
%correct event coding function. event coding function should output a
%structure with  nx1 fields with the names that will be assigned to the nwb
%field with the .task name.

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
    event_data{1} = taskPassiveGLOv1(codes, times);
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