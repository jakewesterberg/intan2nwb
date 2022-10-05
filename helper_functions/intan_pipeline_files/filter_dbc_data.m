tic;
folder_path = 'D:\Intan Recordings\LEO\TEST_DBC128\TEST128DBC_220915_111313';
folder_intan = dir(folder_path);
load('dbc_order.mat')
intan_data = read_intan_header(folder_intan(1).folder);
sampling_rate = intan_data.sampling_rate;
channel_information = intan_data.amplifier_channels;
NUM_CHANNELS = length(intan_data.amplifier_channels);
num_samples = length(intan_data.time_stamp);
downsample_size = length(downsample(intan_data.time_stamp,30));
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

parfor (ii = 1:NUM_CHANNELS,8)
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
spike_data = spike_data(dbc_order+1,:);
low_pass_data = low_pass_data(dbc_order+1,:);
toc;

tic;
%Remove comment to save when satisfactory
save("LEO_test_DBC_data_0909_2022.mat","spike_data","low_pass_data",'-v7.3');
toc;