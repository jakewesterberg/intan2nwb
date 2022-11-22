function nwb = multiNWB(nwb_dirs)

ij_ctr = 1;
for ii = 1 : numel(nwb_dirs)

    dir_info = dir(nwb_dirs{ii});
    for jj = 1 : numel(dir_info)

        if strcmp(dir_info(jj).name(1), '.') || strcmp(dir_info(jj).name(1), 'm')
            continue
        end

        nwb{ij_ctr} = nwbRead([dir_info(jj).folder filesep dir_info(jj).name]);

        ij_ctr = ij_ctr + 1;

    end
end

end