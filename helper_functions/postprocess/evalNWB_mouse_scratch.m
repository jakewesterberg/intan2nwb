num_units = 0;
units_nwb_idx = [];
units_tbl_idx = [];
for ii = 1 : numel(nwb)

    units_per_nwb(ii) = sum(strcmp(nwb{ii}.units.vectordata.get('quality').data(:),'good')); 
    num_units = num_units + units_per_nwb(ii);
    units_nwb_idx = [units_nwb_idx repmat(ii, 1, units_per_nwb(ii))];
    units_tbl_idx = [units_tbl_idx find(strcmp(nwb{ii}.units.vectordata.get('quality').data(:),'good')).'];

end

idx_units = 624;

%%

figure;
%num_channels = numel(nwb{units_nwb_idx(idx_units)}.general_extracellular_ephys_electrodes.vectordata.get('probe_id').data(:));
%nwb_waves = reshape(nwb{units_nwb_idx(idx_units)}.units.waveform_mean.data(:,:), 82, num_channels, []);
%peak_channel = nwb{units_nwb_idx(idx_units)}.units.vectordata.get('peak_channel_id').data(units_tbl_idx(idx_units))+1;
%plot(nwb_waves(:,peak_channel,units_tbl_idx(idx_units)))

u_dat = nwb{units_nwb_idx(idx_units)}.processing.get('convolved_spike_train').nwbdatainterface.get('convolved_spike_train_data').data(units_tbl_idx(idx_units),:);
ts = nwb{units_nwb_idx(idx_units)}.processing.get('convolved_spike_train').nwbdatainterface.get('convolved_spike_train_data').timestamps(:);
evs = nwb{units_nwb_idx(idx_units)}.intervals.get('passive_glo').start_time.data(:);
p45s = nwb{units_nwb_idx(idx_units)}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==45 | ...
    nwb{units_nwb_idx(idx_units)}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==135;

idxs = nearestIdx(ts, evs(p45s));
m_dat = pullVecs(u_dat, idxs, [50 550]); 
plot(-50:550, squeeze(baseline_correct(m_dat,1:50)), 'color', 'k') 
hold on;
yyaxis right
plot(-50:550, squeeze(mean(baseline_correct(m_dat,1:50),3)), 'color', 'r', 'linewidth', 2)
set(gca, 'xlim', [-50 550])