% Question A.) 

EEG = load('data_sets/HW1_Materials/sampleEEGdata.mat').EEG;

elocFile = 'data_sets/HW1_Materials/eloc64C2.txt';

% Extract epochs from 0 to 800 ms. Compute the event-related potentials (ERPs) at each
% electrode. Get the average over all the trials. Select nine time points at which to show
% topographical plots (e.g., 0 to 800 ms in 100ms steps). In one figure, make a series of
% topographical plots at these time points. To increase the signal-to-noise ratio, make each plot
% show the average of activity from 20 ms before until 20 ms after each time point. For example,
% the topographical plot from 200 ms should show average activity from 180 ms until 220 ms.
% Indicate the center time point in a title on each subplot.

epochWindow = [-22 890]; 

% Convert time to sample indices
% We are using deltaT * frequency to get the number of samples
start_index = round((epochWindow(1) - EEG.times(1)) / 1000 * EEG.srate) + 1;
end_index = round((epochWindow(2) - EEG.times(1)) / 1000 * EEG.srate) + 1; 

% epoching our EEG matrix
epochedEegMatrix = EEG.data(:, start_index:end_index, :); 
epochedTimeMatrix = EEG.times(start_index:end_index); 

% Compute the mean over the time dimension (dimension 3)
% I assume this is practically applying a lowpass filter on the data.
% Since all the noise is random, we can assume noise recorded is 
% different in phase. This means when averaged over multiple trials 
% it is suppressed. Therefor the final 64x(len(data))x1 contains 
% the true signal. 
erp = mean(epochedEegMatrix, 3);

%% PLOTTING

% Defining the time points we want to take
start_value = 0; 
step_size = 100; 
num_elements = 9; 
time_points = start_value + (0:(num_elements-1)) * step_size;

% Define the layout for the subplots (adjust rows and columns as needed)
num_rows = 3;
num_cols = ceil(length(time_points) / num_rows);

% In order to show the average activity from 20ms before and 
% 20 ms after, we need to first compute how many points 20ms is
samplesFor20ms = round((20 / 1000) * EEG.srate) + 1;

% Define the desired figure width and height in inches
figureWidth = 15;  % Adjust as needed
figureHeight = 6;  % Adjust as needed

% Create a new figure with the specified size
figure('Units', 'inches', 'Position', [0, 0, figureWidth, figureHeight]);

for i = 1:length(time_points)
    % Find the index corresponding to the current time point
    time_point = time_points(i);
    time_point_index = round(((time_point - epochedTimeMatrix(1)) / 1000) * EEG.srate) + 1; 
    
    % Finding the timepoints index we must average around
    times = [time_point_index - samplesFor20ms , time_point_index + samplesFor20ms]; 

    % Extract the first 4 elements of times of each electrode and calculate their average
    averages = mean(erp(:, times(1):times(2)), 2);
%     averages = erp(:,time_point_index); 
   
    % Create a subplot
    subplot(num_rows, num_cols, i);
    
    % Plot the topographic map
    topoplot(double(averages),elocFile,'eeg'); 
    title([num2str(time_point) 'ms']);
  
    c = colorbar;
    c.Label.String = 'Amplitude in mv';
    caxis([-10 10]);
    
%     sgtitle('ERP data for averaged electrodes +- 20ms') 
end
set(findall(gcf, 'Type', 'axes'), 'FontSize', 12); % Adjust the font size as needed
sgtitle('ERP data averaged +- 20ms, plotted by a center time point');
colormap('Jet')


