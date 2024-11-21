windowSize = 1 * EEG.srate; % 1 second window size
nfft = windowSize; 
overlap = windowSize / 2; % 50 percent overlap

[psd, ~]  = pwelch(win_mean(1,:), windowSize, overlap, nfft, EEG.srate);
