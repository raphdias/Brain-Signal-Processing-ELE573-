params.os.dirs.files = '/Users/raph/Desktop/ele537/matlab/data_sets';
load([params.os.dirs.files '/sampleEEGdata.mat'])

channelData = squeeze(EEG.data(48,:,:))';

% create wavelet
time = -1:1/EEG.srate:1;
f = 5; % frequency of sine wave in Hz
sine_wave = exp(1i*2*pi*f.*time);
s = 4/(2*pi*f); 
gaussian_win = exp(-time.^2./(2*s^2));
wavelet = sine_wave .* gaussian_win;
% half of the wavelet size, useful for chopping off edges after convolution.
halfwaveletsize = ceil(length(wavelet)/2);

% convolve with data
% compute Gaussian
n_conv = length(wavelet) + EEG.pnts - 1;


for itrial = 1:EEG.trials
    eegdata = channelData(itrial,:);
    fft_w = fft(wavelet,n_conv);
    fft_e = fft(eegdata,n_conv);
    ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
    wavelet_conv_data(itrial,:) = ift(halfwaveletsize:end-halfwaveletsize+1);
end

wavelet_conv_erp = mean(wavelet_conv_data,1); 


% create filter and apply to data 
% (more on how to interpret this code in a few chapters!)
nyquist       = EEG.srate/2;
transition_width = 0.2; % percent
filter_low    = 2; % Hz
filter_high   = 8; % Hz
ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(round(1*(EEG.srate/filter_low)),ffrequencies,idealresponse);
for itrial = 1:EEG.trials
    eeg_4to8(itrial,:) = filtfilt(filterweights,1,double(channelData(itrial,:)));
end

eeg_filter_erp = mean(eeg_4to8,1); 

% now plot all the pieces
figure

plot(EEG.times,mean(channelData,1))
hold on
plot(EEG.times,wavelet_conv_erp,'r','linew',2)
plot(EEG.times,eeg_filter_erp,'m','linew',2)
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
legend({'Raw data';'wavelet convolved';'band-pass filtered'})
