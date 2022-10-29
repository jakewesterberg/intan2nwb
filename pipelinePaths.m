function pp = pipelinePaths(varargin)

% Defaults
pp.RAW_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_0_RAW_DATA\';
pp.CAT_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_1_CAT_DATA\';
pp.BIN_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_2_BIN_DATA\';
pp.SPK_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_3_SPK_DATA\';
pp.SSC_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_4_SSC_DATA\';
pp.CNX_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_5_CNX_DATA\';
pp.NWB_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_6_NWB_DATA\';
pp.EPO_DATA     = 'Z:\_DATA\_BL_DATA_PIPELINE\_7_EPO_DATA\';

pp.CONDA        = 'C:\Users\westerja.VUDS\Anaconda3';
pp.REPO         = 'C:\Users\westerja.VUDS\Documents\GitHub\intan2nwb\';

pp.SCRATCH      = [userpath filesep 'scratch'];

% Varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'RAW_DATA'}
            pp.RAW_DATA = varargin{varStrInd(iv)+1};
        case {'CAT_DATA'}
            pp.CAT_DATA = varargin{varStrInd(iv)+1};
        case {'BIN_DATA'}
            pp.BIN_DATA = varargin{varStrInd(iv)+1};
        case {'SPK_DATA'}
            pp.SPK_DATA = varargin{varStrInd(iv)+1};
        case {'SSC_DATA'}
            pp.SSC_DATA = varargin{varStrInd(iv)+1};
        case {'CNX_DATA'}
            pp.CNX_DATA = varargin{varStrInd(iv)+1};
        case {'NWB_DATA'}
            pp.NWB_DATA = varargin{varStrInd(iv)+1};
        case {'EPO_DATA'}
            pp.EPO_DATA = varargin{varStrInd(iv)+1};
        case {'CONDA'}
            pp.CONDA = varargin{varStrInd(iv)+1};
        case {'REPO'}
            pp.REPO = varargin{varStrInd(iv)+1};
        case {'SCRATCH'}
            pp.SCRATCH = varargin{varStrInd(iv)+1};
    end
end

end