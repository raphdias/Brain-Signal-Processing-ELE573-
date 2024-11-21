% Raphael Dias Homework 2

%% Question A.) 

% Load Participant Data
clear; 
clc; 

% Define the folder path and load eloc data into a table
folderPath = 'data_sets';
eloc16C2 = readtable(fullfile(folderPath,'eloc16C2.txt')); 

% Create a cell array to store EEG data for each subject
% Use the 'dir' function to list all files in the directory
fileList = dir(fullfile(folderPath, '*.mat')); 
numSubjects = numel(fileList); % number of subjects
eegDataCell = cell(1, numSubjects);

% [subject, open, closed, f, peakFreqAtAlpha] 
psdDataCell = cell(numSubjects,5); 

% Load or generate EEG data for each subject
% For this example, let's generate random data
for subject = 1:numSubjects
    % Generate random EEG data (replace with your data loading code)
    eegData = load(fullfile(folderPath, fileList(subject).name)).data;
    
    % Store EEG data for the subject
    eegDataCell{subject} = eegData;
end 


%% Computing welch 
% The idea of welches method - Isolate
% a little piece of the data, get the fourier transform and compute the 
% power spectrum of the snippet
%
% Welchs method has an overlap, this is where the next window starts in 
% comparison to where you first signal ends. This is done because 
% of a tapering that happens, and is attenuated when overlapped
%
% in the end, we average the power spectra of each snippet together
% Welchs method genreally smooths out the power spectrum 

% Define the parameters
OzIndex = 16; % Oz is stored in column 16 
Fs = 256; % hz
windowSize = 256; % 1 second window size
nfft = windowSize; 
overlap = windowSize / 2; % 50 percent overlap 

figure(1); 
for subject = 1:numSubjects
    
    % defining eyes open and closed data for channel Oz
    eyesOpenData   = eegDataCell{subject}.EyesOpen(:,16); 
    eyesClosedData = eegDataCell{subject}.EyesClosed(:,16); 

    [psdOpen, f] = pwelch(eyesOpenData, windowSize, overlap, nfft, Fs);

    [psdClosed, f_] = pwelch(eyesClosedData, windowSize, overlap, nfft, Fs);

    psdDataCell{subject,1} = fileList(subject).name;
    psdDataCell{subject,2} = eyesOpenData; 
    psdDataCell{subject,3} = eyesClosedData; 
    psdDataCell{subject,4} = f; 

    % identifying the peak frequency occurance at alpha (8 - 12hz); 
    [M_, IClosed] = max(psdClosed(8:12)); 
    psdDataCell{subject,5} = IClosed + 7; 

    % Plot PSD in a subplot
    subplot(3, 2, subject);

    % Plot PSD in a subplot
    semilogy(f, psdOpen, 'b', 'LineWidth', 1.5); % Eyes open condition in blue
    hold on;
    semilogy(f, psdClosed, 'r', 'LineWidth', 1.5); % Eyes closed condition in red
    xlim([0 70]);
    ylim([0 100]);

    % Add labels and legend
    title(strrep(fileList(subject).name,'_','-'))
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
    legend('Eyes Open', 'Eyes Closed');
end

sgtitle('PSD using 256-pt pwelch with 50% overlap')



%% Question B.)  (BROKEN)

% init our vector for topoplot
% Define the variable names and data types
variableNames = {'SubjectName', 'EyesOpen','EyesClosed'};
variableTypes = {'string', 'single','single'};

% Create an empty table with variable names and data types
subjectFrequencyTable = table('Size', [0, numel(variableNames)], 'VariableNames', variableNames, 'VariableTypes', variableTypes);

% Here we are going to pwelch each channel, and then 
% gather the amplitude at our maximum frquency: pdDataCell{subject,5} 
f = figure;
tcl = tiledlayout(f,"flow",TileSpacing="tight");

elocFile = 'data_sets/eloc16C2.txt'; 
for subject = 1:numSubjects

    % init vectors
    openAmplitudeVector = []; 
    closedAmplitudeVector = []; 

    for channel = 1:height(eloc16C2)

        % data from channel given a condition
        eyesOpenData   = eegDataCell{subject}.EyesOpen(:,channel); 
        eyesClosedData = eegDataCell{subject}.EyesClosed(:,channel);
        
        % psd of each condition
        [psdOpen, f] = pwelch(eyesOpenData, windowSize, overlap, nfft, Fs);
        [psdClosed, f_] = pwelch(eyesClosedData, windowSize, overlap, nfft, Fs);

        openAmplitudeVector(end+1) = psdOpen(psdDataCell{subject,5});
        closedAmplitudeVector(end+1) = psdClosed(psdDataCell{subject,5});
    end

    % Create a new row to append
    subjectFrequencyTable = [subjectFrequencyTable;{string(fileList(subject).name), openAmplitudeVector,closedAmplitudeVector}];
    
    % plotting here
    nexttile; 
    topoplot(double(subjectFrequencyTable.EyesOpen(subject,:)),elocFile,'eeg')
    title(strrep(fileList(subject).name,'_','-')); 
    nexttile; 
    topoplot(double(subjectFrequencyTable.EyesClosed(subject,:)),elocFile,'eeg')
    title(['Oz frquency - ' num2str(psdDataCell{subject,5}) ' Hz']); 
    cb = colorbar; 
  
end
colormap('jet')

%% Question C.) 

T7 = 3; 
T8 = 4; 

% Define the parameters
OzIndex = 16; % Oz is stored in column 16 
Fs = 256; % hz
windowSize = 256; % 1 second window size
nfft = windowSize; 
overlap = windowSize / 2; % 50 percent overlap 

figure; 

for subject = 1:numSubjects

    % defining eyes open and closed data for channel T7 and T8 
    clenchT7   = eegDataCell{subject}.Clench(:,T7); 
    clenchT8 = eegDataCell{subject}.Clench(:,T8); 

    % psd for channels 
    [psdT7, f] = pwelch(eyesOpenData, windowSize, overlap, nfft, Fs);
    [psdT8, f_] = pwelch(eyesClosedData, windowSize, overlap, nfft, Fs);

    % Plot PSD in a subplot
    subplot(3, 2, subject);

    % Plot PSD in a subplot
    semilogy(f, psdT7, 'b', 'LineWidth', 1.5); % T7 in blue
    hold on;
    semilogy(f, psdT8, 'r', 'LineWidth', 1.5); % T8 condition in red
    xlim([0 70]);
    ylim([0 100]);

    % Add labels and legend
    title(strrep(fileList(subject).name,'_','-'))
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
    legend('T7', 'T8');
end

sgtitle('PSD using 256-pt pwelch with 50% overlap')



%% Question D.) (Do)

frequencies = [10 25 40 65]; 

% Here we are going to pwelch each channel, and then 
% gather the amplitude at our maximum frquency: pdDataCell{subject,5} 
% f = figure;
% tcl = tiledlayout(f,"flow",TileSpacing="tight");

elocFile = 'data_sets/eloc16C2.txt'; 
f = figure;
tcl = tiledlayout(f,"flow",TileSpacing="loose");

for subject = 1:numSubjects

    % init vectors
    clenchVector = []; 

    for channel = 1:height(eloc16C2)

        % data from channel given a condition
        teethClenchedData = eegDataCell{subject}.Clench(:,channel); 
        
        % psd of each condition
        [psdClench, f] = pwelch(teethClenchedData, windowSize, overlap, nfft, Fs);
     
        clenchVector = [clenchVector psdClench(frequencies)];
      
    end

    % Create a new row to append    
    % plotting here
    
    
    for plotIndex = 1:length(frequencies)
        nexttile; 
        topoplot(double(clenchVector(plotIndex,:)),elocFile,'eeg')
        title([num2str(frequencies(plotIndex)) ' Hz']); 
    end
    cb = colorbar; 
end
colormap('jet')


%% Question 2.) 

Fz = 2;

eeg_data_unfiltered_cell = cell(1, numSubjects);
eog_data = cell(1, numSubjects);
bp_filt = cell(1, numSubjects);

for subject=1:numSubjects

    % saving eog data and filtering 1-15hz
    eog_data_unfiltered = double(reshape(squeeze(eegDataCell{1, subject}.Blink(:,Fz)), 1, []));
    eeg_data_unfiltered_cell{subject} = eog_data_unfiltered; 

    % Apply a bandpass filter to focus on frequencies between 1 Hz and 8 Hz
    % this is a butterworth filter with order = 4.
    % We are also splitting the filter into 2 low order filter - applying
    % sequentially.
    filter_order = 4; 
    cutoff_frequency_1 = 1; 
    cutoff_frequency_2 = 8;


    bpfilt = designfilt('bandpassfir', ...
        'FilterOrder',filter_order,'CutoffFrequency1',cutoff_frequency_1, ...
        'CutoffFrequency2',cutoff_frequency_2,'SampleRate',Fs);

    % zero phase filtered on the data
    eog_filtered = filtfilt(bpfilt, eog_data_unfiltered); 
    % saving just to plot
    bp_filt{subject} = eog_filtered; 

    % Detrend the absolute value signal, removing the constant (DC) component
    eog_data{subject} = detrend(eog_filtered, 'constant').^2;
end

%% Plot filtered singal compared to orignal signal

% example of each step 
subject = 1; 
figure; 
subplot(3,1,1);
plot(eeg_data_unfiltered_cell{subject});
title("Oringal Signal")
subplot(3,1,2);
plot(bp_filt{subject}); 
title("Bandpass Filter 1-8hz")
subplot(3,1,3)
plot(eog_data{subject});
title("Detrend^2: Removed constant DC component and sqaured the signal")


%% detecting eyebliks and saving indexes to array 

spikesCell = cell(1,numSubjects); 

for subject=1:numSubjects

    eog_filt = eog_data{subject};

    % threshold is the mean plus the standard deviation
    threshold = mean(eog_filt) + std(eog_filt); 

    eblength = round(Fs*0.25); %length of eyeblink(200 ms) in samples;
    spikes = [];

    for i = eblength:length(eog_filt)-eblength
    if abs(eog_filt(i))>threshold && ... %bigger than threshold
       all(abs(eog_filt(i))>=abs(eog_filt(i-eblength+1:i+eblength))) %biggest in surrounding 400ms
        spikes = [spikes i];
    end

    spikesCell{subject} = spikes; 
    end         
end

fprintf('Eyeblink Spikes Detected %d\n', length(spikes));

%% Plotting original signal with eyeblinks marked

startPoint = 1; 
endPoint = 8000; 
% spikeLessThanCondition = (spikes >= startPoint & spikes <= endPoint);
% spikeIndexes = find(spikeLessThanCondition); 

figure; 
for i=1:numSubjects

    subplot(5,1,i); 
    plot(eeg_data_unfiltered_cell{i});
    hold on; 
    xline(spikesCell{i},'r--'); % 'r--' specifies red dashed lines
    title(strrep(fileList(i).name,'_','-'))
end
sgtitle('Eye Blink Detection Algorithim using Blink Condition - Accross all subjects')



