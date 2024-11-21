% ProbB
% Generate 5 sine waves b/w 2-30Hz
close all;
clear;
Fs = 256;  % samples per second
dt = 1/Fs;   % seconds per sample (time increment)
t = (-1:dt:1); % 5 seconds of data (time)
t = t(1:Fs*2);
sine_waves = zeros(5,512);%sin(2*pi*f*t); 

load('sampleEEGdata.mat');
elec_file_path = 'eloc64C2.txt';

% Open the file for reading
fileID = fopen(elec_file_path, 'r');

% Read the data using textscan
data_elec_file = textscan(fileID, '%d %f %f %s', 'Delimiter', '\t');


% Close the file
fclose(fileID);

% Extract the electrode number, coordinates, and name
electrode_number = data_elec_file{1};
coordinates_x = data_elec_file{2};
coordinates_y = data_elec_file{3};
electrode_name = data_elec_file{4};
ch_count = EEG.nbchan;
num_trials = EEG.trials;
time_points = EEG.times;
data_eeg = EEG.data(:,:,1);

% Frequencies of interest
frequencies = 2:7:30;

% Create and plot the sine waves
figure;

for i = 1:length(frequencies)
    f = frequencies(i);
    y = sin(2*pi*f*t);
    sine_waves(i,:) = y(1:Fs*2);
    
    subplot(length(frequencies), 1, i);
    plot(t, y);
    title(['Sine wave of ' num2str(f) ' Hz']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    ylim([-1, 1]);
end

%Create gaussian window
a = 1; % Amplitude
m = 0; % Offset
n = 4; % Num of cycles
GaussWins = zeros(5,512);

figure;
for i = 1:length(frequencies)
    f = frequencies(i);
    s = n/(2*pi*f);
    y = a.*(exp(-(t-m).^2/(2*s^2)).*exp(j*2*pi*f*t));
    GaussWins(i,:) = y(1:Fs*2);
    subplot(length(frequencies), 1, i);
    plot(t, GaussWins(i,:));
    title(['Gauss window of ' num2str(f) ' Hz']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    ylim([-1, 1]);
end


[num_channels, num_timepoints, trials] = size(data_eeg);
[num_wavelets, wavelet_length] = size(GaussWins);
num_timepoints = 512; 

% Initialize an empty matrix to store the convolved data
convolved_data = zeros(num_channels, num_timepoints, num_wavelets);

% Amount to be truncated from each end of the convolved data
truncate_length = floor(wavelet_length / 2);

for i = 1:num_wavelets
    for k = 1:num_channels
        % Convolve the EEG data with the current wavelet

        % changing to fft convolution
        freq_conv = ifft(fft(data_eeg(k,1:512)).*fft(GaussWins(i,:)));

%         temp_conv = conv(data_eeg(k, :), GaussWins(i, :), 'same');
        % Truncate the convolved data to match the original data length
        convolved_data(k, :, i) = freq_conv;%temp_conv(truncate_length:end - truncate_length);
    end
end

% Initialize the matrix to store power and phase
% Dimensions: [num_timepoints, num_wavelets, num_channels, 2]
power_phase_matrix = zeros(num_timepoints, num_wavelets, num_channels, 2);

for i = 1:num_wavelets
    for k = 1:num_channels
        % Extract complex data
        complex_data = convolved_data(k, :, i);
        
        % Compute the power (magnitude squared)
        power = abs(complex_data).^2;
%         disp(min(power))
%         disp(max(power))
        
        % Compute the phase (angle)
        phase = angle(complex_data);
        
        % Store in the matrix
        power_phase_matrix(:, i, k, 1) = power;
        power_phase_matrix(:, i, k, 2) = phase;
    end
end
time_event = 180;
time_window_start = time_event - 20;
time_window_end = time_event + 20;
ind_time_window = find(time_points >= time_window_start & time_points <= time_window_end);
%find(time_points >= 0.12 & time_points <= 0.180);
% n=0.36;
% [val,ind_time_window]=min(abs(time_points-n));

% Calculate the average activity within the time window for all channels
power_plot = double(squeeze(mean(power_phase_matrix(ind_time_window,:,: , 1),1)));
phase_plot = double(squeeze(mean(power_phase_matrix(ind_time_window,:,: ,2),1)));

% Calculate the global minimum and maximum for power and phase
power_max = max(power_plot(:));
power_min = min(power_plot(:));

% Calculate the global minimum and maximum for power 
power_log_min = log(power_min + 1); 
power_log_max = log(power_max + 1);
figure;
% Iterate through the 5 wavelet frequencies and plot power and phase
for freq = 1:5
    % Power plot
    subplot(1, 5, freq);
    power_plot_freq = double(squeeze(mean(power_phase_matrix(ind_time_window,freq,:,1),1)));
    topoplot(log(power_plot_freq + 1), 'eloc64C2.txt', 'eeg'); % Taking log for exponential scaling
    title(['Frequency: ' num2str(frequencies(freq)) ' Hz (Power)']);
    caxis([power_log_min, power_log_max]);
end

% Create exponential colormap for power
num_colors = 256; % Number of colors in the colormap
exp_colormap = jet(num_colors);
colormap(exp_colormap);

% Adjusting the colorbar for exponential scale
cbar = colorbar('Position', [0.93, 0.1, 0.02, 0.8]);
% Define ticks for colorbar in log scale
cbar_ticks = linspace(power_log_min, power_log_max, 5);
set(cbar, 'Ticks', cbar_ticks);
% Convert log scale ticks back to original scale for labeling
cbar_ticklabels = round(exp(cbar_ticks) - 1);
set(cbar, 'TickLabels', cbar_ticklabels);

set(findall(gcf, 'Type', 'axes'), 'FontSize', 12); % Adjust the font size as needed


phase_max =  max(phase_plot(:));
phase_min = min(phase_plot(:));
% Create a figure to hold the subplots for phase
figure;
title('Topographical Plots of Phase for Different Frequencies');

% Iterate through the 5 wavelet frequencies and plot power and phase
for freq = 1:5
    
    % Phase plotc
    subplot(1, 5, freq); % 2x5 grid, bottom row
    phase_plot = double(squeeze(mean(power_phase_matrix(ind_time_window,freq,:,2),1)));
    topoplot(phase_plot, 'eloc64C2.txt', 'eeg');
    title(['Frequency: ' num2str(frequencies(freq)) ' Hz (Phase)']);
    caxis([phase_min, phase_max]);
end
colormap('jet');
colorbar('Position', [0.93, 0.1, 0.02, 0.8]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 12); % Adjust the font size as needed
