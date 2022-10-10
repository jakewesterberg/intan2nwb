function val = returnGSNum(dat, ii, jj)
if strcmp(class(dat), 'double')
    val = dat(ii);
elseif strcmp(class(dat), 'cell')
    temp_array_1 = strtrim(split(dat{ii}, ','));
    val = str2double(temp_array_1{jj});
    clear temp_array_1
end
end