% Loading in the data set
params.os.dirs.files = '//Users/raph/Desktop/Fall 2024/ele537/matlab/data_sets';
load([params.os.dirs.files '/sampleEEGdata.mat'])

t = -1:1/EEG.srate:1; 
freq_vector=[4,30];
Number_of_freq=length(freq_vector);
n_cycles=5;
% you can also use linspace as follows:
% freq_vector=round(linspace(min_freq, max_freq, 5));
% Then you will loop on each frequency in order to create your complex sine
% wave and Gaussian and multiply them together and store the result in a
% matrix "Morlet"


% In order to do MORLET Wavelet and not Complex Morlet Wavelet 
% This Loop Is What You Want To Change
for i=1:Number_of_freq
    f=freq_vector(i);
    s=n_cycles/(2*pi*f);
    A=1/sqrt(s*pi^0.5);%complex sine wave amplitude
    ComplexSineWave=A*exp(1i*2*pi*f*t);
    %Gaussian wave , 4 cycles
    GW=exp(-1*t.^2/(2*s.^2));
    Morlet(i,:)=ComplexSineWave.*GW;%creating our wavelet by element wise
end

% extract data related to channel Cz
ch = 48;
channelData=squeeze(EEG.data(ch,:,:))';

% half of the wavelet size, useful for chopping off edges after convolution.
halfwaveletsize = ceil(length(Morlet(1,:))/2);

% convolve with data
% compute Gaussian
n_conv = length(Morlet(1,:)) + EEG.pnts - 1;
for i = 1:Number_of_freq
    for itrial = 1:EEG.trials
        eegdata = channelData(itrial,:);
        fft_w = fft(Morlet(i,:),n_conv);
        fft_e = fft(eegdata,n_conv);
        ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
        wavelet_conv_data(i,itrial,:) = ift(halfwaveletsize:end-halfwaveletsize+1);
    end
end

%%

wavelet_conv_erp = squeeze(mean(wavelet_conv_data,2)); 

clear i itral 

Power=abs(wavelet_conv_erp).^2;
Real=real(wavelet_conv_erp);


%%

figure;
subplot(221)
plot(EEG.times,Power(1,:)); 
ylabel('Power')
xlabel('time (ms)')
title('Power at 5hz')
subplot(222)
plot(EEG.times,Real(1,:)); 
ylabel('Real Value')
xlabel('time (ms)')
title('Real at 5hz')
subplot(223)
plot(EEG.times,Power(2,:)); 
ylabel('Power')
xlabel('time (ms)')
title('Power at 25hz')
subplot(224)
plot(EEG.times,Real(2,:)); 
ylabel('Real Value')
xlabel('time (ms)')
title('Real at 25hz')

sgtitle('CMW Convolution')


