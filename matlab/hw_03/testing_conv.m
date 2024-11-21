clear; 
close all; 

load('sampleEEGdata.mat')
t = -1:1/256:1; 

% two signals of same length
eeg_data = EEG.data(1,1:513,1);
wav = morlet(t,10); 

% using conv function 
timeConv = conv(eeg_data, wav, 'same'); 

% using fft
eeg_fft = fft(eeg_data);
wav_fft = fft(wav); 
freqConv = ifft(eeg_fft.*wav_fft); 

figure(); 
plot(log(abs(timeConv).^2), '-o'); 
hold on; 
plot(log(abs(freqConv).^2), '--'); 

clear t






function wavelet = morlet(t, f)
    cycles = 4; 
    sigma = cycles / (2 * pi * f);
    
    % To make a Morlet wavelet create a sine wave, create a Gaussian wave, and
    % multiply them together.
    freqBandScale = 1/(sigma*sqrt(pi))^(1/2); 
    wavelet = freqBandScale.*exp(1i * 2 * pi * f * t) .* exp(-t.^2 / (2 * sigma^2));
end

