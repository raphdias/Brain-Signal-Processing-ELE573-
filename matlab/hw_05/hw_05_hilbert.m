clear; close all; 

% Loading in the data set
params.os.dirs.files = '/Users/raph/Desktop/ele537/matlab/data_sets';
load([params.os.dirs.files '/sampleEEGdata.mat'])
elocFile = [params.os.dirs.files '/HW1_Materials/eloc64C2.txt'];

%%This script guides you through filtering the data prior to implementing the hilbert. Filter weights are designed using FIR1 filter and are given to matlab filtfilt function. 
% For more info please read matlab fir1 and filtfilt functions
%
%
%
freqs=[5 25]; %center frequency; 5 or 25 Hz
width=3; %frequency band width
ch = 48;%put the channel number corresponding to channel "Cz". You need to take a look at the eloc64.txt file to find the corresponding channel number.
for ifreq = 1:length(freqs)

    freq = freqs(ifreq); 

    Freqmin=freq-width; %miniumum frequency
    Freqmax=freq+width; %maximum frequency
    % Note - I changed filter order to fit the bw of 3hz
    Filterorder=round(1*(EEG.srate/Freqmin)); %length of the filter in time domain at least 3 times of the lowest activity (lowest frequency in Hz)
    FWeights= fir1(Filterorder,[Freqmin/(EEG.srate/2) Freqmax/(EEG.srate/2)], 'bandpass');
    ChannelData=squeeze(EEG.data(ch,:,:));
    ChannelData=ChannelData';%trials*samples
    for trial=1:EEG.trials %number of trials
        FilteredData(trial,:)=filtfilt(FWeights,1,double(ChannelData(trial,:)));
    end

    % for each trial compute the fft
    for itrial = 1:EEG.trials
        f(itrial,:) = fft(FilteredData(itrial,:));
    end

    % Create complex copy of f
    complexf = 1i*f; 
    % Idetifying positive and negative frequencies
    % positive frequencies (pf) = 0<pf<nyquist
    % negative frequencies (nf) = nyquist<pf
    
    % Amplitudes in real valued signals get split 
    % between positive frequencies and negative frequencies. 
    % The frequencies are split evenly in the fourier transform, meaning
    % the first half of the fft is positive, and second half is negative 
    posF = 2:floor(EEG.pnts/2)+mod(EEG.pnts,2);
    negF = ceil(EEG.pnts/2)+1+~mod(EEG.pnts,2):EEG.pnts;
    
    % rotate -90 for positive and +90 for negative, then add to f
    f(:,posF) = f(:,posF) + -1i*complexf(:,posF);
    f(:,negF) = f(:,negF) +  1i*complexf(:,posF);
    
    hilbert_manual(ifreq,:,:) = mean(ifft(f,[],2),1)';
    % compare with Matlab function hilbert
    hilbertm(ifreq,:,:) = hilbert(mean(FilteredData,1)');
end
% % 
% plot results
figure; 
count = 1; 
for i = 1:length(freqs)
    subplot(2,2,count)
    plot(EEG.times,abs(hilbertm(i,:)).^2)
    ylabel('Power')
    xlabel('time (ms)')
    hold on
    plot(EEG.times,abs(hilbert_manual(i,:)).^2,'r')
    ylabel('Power')
    xlabel('time (ms)')
    legend({'Matlab Hilbert function';'"manual" Hilbert'})
    title(['Power of Hilbert transform at ' num2str(freqs(i)) ' Hz'], 'Interpreter','none')
    count = count + 1; 
    
    subplot(2,2,count)
    plot(EEG.times,real(hilbertm(i,:)))
    ylabel('Real Value')
    xlabel('time (ms)')
    hold on
    plot(EEG.times,real(hilbert_manual(i,:)),'r')
    ylabel('Real Value')
    xlabel('time (ms)')
    legend({'Matlab Hilbert function';'"manual" Hilbert'})
    title(['Real of Hilbert transform at ' num2str(freqs(i)) ' Hz'], 'Interpreter','none')
    count = count + 1; 

end

sgtitle('Hilbert Transform')

