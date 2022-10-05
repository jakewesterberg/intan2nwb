tic;
folder_path = 'C:\Users\Recording\Documents\Intan Recording\CurrentRecording\LEO_GLO_CONTROL\GLO_CONTROL_VPROBE_220908_220908_113515';
folder_intan = dir(folder_path);
intan_data = read_intan_header(folder_intan(1).folder);
load('v_probe_order.mat')
sampling_rate = intan_data.sampling_rate;
channel_information = intan_data.amplifier_channels;
NUM_CHANNELS = length(intan_data.amplifier_channels);
num_samples = length(intan_data.board_adc_data(1,:));
downsample_size = length(downsample(intan_data.board_adc_data(1,:),30));
spike_data = zeros(NUM_CHANNELS,downsample_size);
low_pass_data = zeros(NUM_CHANNELS,downsample_size);
toc;

tic;
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',500,'CutoffFrequency2',5000, ...
         'SampleRate',sampling_rate);
lpFilt = designfilt('lowpassfir','FilterOrder',20, ...
         'CutoffFrequency',250, ...
         'SampleRate',sampling_rate);
parfor (ii = 1:NUM_CHANNELS,4)
    current_fid = fopen(folder_path+"\amp-" + channel_information(ii).native_channel_name + ".dat");
    current_data = gpuArray(double(fread(current_fid,num_samples,'int16')) * 0.195);
    fclose(current_fid);
    spike_data(ii,:) = gather(downsample(filtfilt(lpFilt,abs(filtfilt(bpFilt,current_data))),30));
    low_pass_data(ii,:) = gather(downsample(filtfilt(lpFilt,current_data),30));
end
toc;

%Rearrange the channels to the order on the probe (starts at 0, +1 so it
%matches matlab indexing)
tic;
spike_data = spike_data(v_probe_order,:);
low_pass_data = low_pass_data(v_probe_order,:);
toc;

tic;
%Remove comment to save when satisfactory
save("LEO_glo_v_probe_data_0830_2022_third_test.mat","spike_data","low_pass_data",'-v7.3');
toc;