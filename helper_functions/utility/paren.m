function dat = paren(dat, varargin)
if isstr(varargin{1})
    if strcmp(varargin{1}, 'vec')
        dat = dat(:);
    end
else
    dat = dat(varargin{:});
end
end