num_units = 0;
units_nwb_idx = [];
units_tbl_idx = [];
for ii = 1 : numel(nwb)

    units_per_nwb(ii) = sum(strcmp(nwb{ii}.units.vectordata.get('quality').data(:),'good')); 
    num_units = num_units + units_per_nwb(ii);
    units_nwb_idx = [units_nwb_idx repmat(ii, 1, units_per_nwb(ii))];
    units_tbl_idx = [units_tbl_idx find(strcmp(nwb{ii}.units.vectordata.get('quality').data(:),'good')).'];

end

lfp_nwb_idx = [];
lfp_tbl_idx = [];
for ii = 1 : numel(nwb)

    lfp_per_nwb(ii) = size(nwb{ii}.acquisition.get('probe_0_lfp').electricalseries.get('probe_0_lfp_data').data(:,:),1); 
    lfp_nwb_idx = [lfp_nwb_idx repmat(ii, 1, lfp_per_nwb(ii))];
    lfp_tbl_idx = [lfp_tbl_idx 1:lfp_per_nwb(ii)];

end


u_dat = nwb{idx_session}.acquisition.get('photodiode_1_tracking').timeseries.get('photodiode_1_tracking_data').data(:);
[bwb, bwa] = butter(2, 75/625, 'low');
u_dat = filtfilt(bwb, bwa, u_dat);

figure;
ts = nwb{idx_session}.acquisition.get('photodiode_1_tracking').timeseries.get('photodiode_1_tracking_data').timestamps(:);
evs = nwb{idx_session}.intervals.get('passive_glo').start_time.data(:);
pAA = nwb{idx_session}.intervals.get('passive_glo').vectordata.get('correct').data(:)==1 & ...
    (nwb{idx_session}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==45 | ...
    nwb{idx_session}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==135);
idxs = nearestIdx(ts, evs(pAA));
m_dat = pullVecs(u_dat', idxs, [50 550]); 
plot(-50:550, squeeze(baseline_correct(m_dat,1:50)), 'color', 'k') 
hold on;
yyaxis right
plot(-50:550, squeeze(mean(baseline_correct(m_dat,1:50),3)), 'color', 'r', 'linewidth', 2)
set(gca, 'xlim', [-10 150]) 


idx_session = 4;
probe_no = 'probe_0_muae';

%u_dat = common_average_reference(nwb{idx_session}.acquisition.get(probe_no).electricalseries.get([probe_no '_data']).data(:,:));
u_dat = nwb{idx_session}.acquisition.get(probe_no).electricalseries.get([probe_no '_data']).data(:,:);
for i = 1:size(u_dat,1)
    u_dat(i,:) = smooth(u_dat(i,:),50);
end
u_chan = size(u_dat,1);

idx_lfp = 27;
figure;

subplot(1,3,1)
ts = nwb{idx_session}.acquisition.get(probe_no).electricalseries.get([probe_no '_data']).timestamps(:);
evs = nwb{idx_session}.intervals.get('passive_glo').start_time.data(:);
% pAA = nwb{idx_session}.intervals.get('passive_glo').vectordata.get('correct').data(:)==1 & ...
%      (nwb{idx_session}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==45 | ...
%      nwb{idx_session}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==135);
% 
% p45 = nwb{idx_session}.intervals.get('passive_glo').vectordata.get('correct').data(:)==1 & ...
%      nwb{idx_session}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==45;
% p135 = nwb{idx_session}.intervals.get('passive_glo').vectordata.get('correct').data(:)==1 & ...
%      nwb{idx_session}.intervals.get('passive_glo').vectordata.get('orientation').data(:)==135;

p45 = nwb{idx_session}.intervals.get('passive_glo').vectordata.get('correct').data(:)==1 & ...
     nwb{idx_session}.intervals.get('passive_glo').vectordata.get('go_gloexp').data(:)==1;
p135 = nwb{idx_session}.intervals.get('passive_glo').vectordata.get('correct').data(:)==1 & ...
     nwb{idx_session}.intervals.get('passive_glo').vectordata.get('gloexp').data(:)==1 & ...
     nwb{idx_session}.intervals.get('passive_glo').vectordata.get('presentation').data(:)==3;

idxs = nearestIdx(ts, evs(p135));
%m_dat = rmOutlier(baseline_correct(pullVecs(u_dat, idxs, [50 550]),1:50));
m_dat = baseline_correct(pullVecs(u_dat, idxs, [50 550]),1:50);
[mean_out, lower, upper] = confidence_interval(m_dat);
hold on;
plot(-50:550,squeeze(m_dat(idx_lfp,:,randi(size(m_dat,3),[1,100]))), 'color', 'k')
plot_ci(lower(idx_lfp,:), upper(idx_lfp,:), -50:550)
plot(-50:550, squeeze(mean_out(idx_lfp,:)), 'color', 'b', 'linewidth', 2)
set(gca, 'xlim', [-50 550])

subplot(1,3,2)
idxs = nearestIdx(ts, evs(p45));
m_dat = baseline_correct(pullVecs(u_dat, idxs, [50 550]),1:50);
[mean_out, lower, upper] = confidence_interval(m_dat);
hold on;
plot_ci(lower(idx_lfp,:), upper(idx_lfp,:), -50:550)
plot(-50:550, squeeze(mean_out(idx_lfp,:)), 'color', 'b', 'linewidth', 2)

idxs = nearestIdx(ts, evs(p135));
m_dat = rmOutlier(baseline_correct(pullVecs(u_dat, idxs, [50 550]),1:50));
[mean_out, lower, upper] = confidence_interval(m_dat);
plot_ci(lower(idx_lfp,:), upper(idx_lfp,:), -50:550)
plot(-50:550, squeeze(mean_out(idx_lfp,:)), 'color', 'r', 'linewidth', 2)
set(gca, 'xlim', [-50 550])

subplot(1,3,3)
imagesc(-50:550, 1:u_chan, mean_out)
set(gca, 'xlim', [-50 550]) 
colorbar