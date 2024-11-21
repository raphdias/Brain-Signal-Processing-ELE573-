params.os.dirs.files = '/Users/raph/Desktop/ele537/matlab/data_sets';
load([params.os.dirs.files '/sampleEEGdata.mat'])


%% Compating to FIR bandpass filter
freqs = [5 25]; 
width = 3;
ch = 48; 

figure; 
for i = 1:length(freqs)
    freq = freqs(i); 
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
    % plot the filter data
    subplot(2,1,i)
    plot(EEG.times,mean(FilteredData,1)'); 
    title(['Bandpass at ' num2str(freq)]); 
    xlabel('time (ms)')
    ylabel('Microvolts (uv)')
end