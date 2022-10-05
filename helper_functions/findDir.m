function dir_cont = findDir(dir_in, search_exp)

dir_str = dir(dir_in);
dir_cont = {};

for itt_str = 1 : length(dir_str)
    
    if dir_str(itt_str).isdir && ~strcmp(dir_str(itt_str).name,'.') && ...
            ~strcmp(dir_str(itt_str).name,'..')
        
        file_name = fullfile(dir_in,dir_str(itt_str).name);
        
        temp_1 = regexp(lower(dir_str(itt_str).name), lower(search_exp), 'once');
        if ~isempty(temp_1)
            dir_cont = [dir_cont; {file_name}];
        end

        temp_dir_cont = findDir(file_name,search_exp);
        
        if ~isempty(temp_dir_cont)            
            dir_cont = [dir_cont; temp_dir_cont];
        end
    end
end
end