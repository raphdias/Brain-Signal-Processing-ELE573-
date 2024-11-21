% Question C.) 

% Repeat step (a) by applying large Laplacian filter. Compare step (c) 
% with step (a) and clearly explain your observations and comments 
% (hint: To obtain the distance and the surrounding electrodes, transfer 
% the polar coordinates in eloc64C2.txt file into Cartesian. Then for each 
% electrode of interest keep those electrodes that are in radius 
% [0.18 0.28], remove the rest, continue obtaining your weights, and then 
% obtain the Laplacian filtered signal).

%% Start


% Background: 
% Laplcian Filter serves as a high-pass filter that enhances localized 
% activity while suppresses the diffusion activity
% '''
% The implementation of the Laplacian in EEG filtering on the voltage at 
% each electrode is to subtract the weighted voltages from the surrounding 
% electrodes from the voltage recording at current electrode, where the 
% weight is electrode distance dependent.
% '''

radius = [0.18 0.28]; 

% Convert eloc64C2.txt file's polar coords into Cartesian
polar_coords = readmatrix('data_sets/HW1_Materials/eloc64C2.txt'); 
cart_coords = []; 

% save the size of polar_coords, should be electrode x coords
n_electrodes = size(polar_coords); 

% converting each polar coordinate to cartesian coords
for electrode_index = 1:n_electrodes(1)

    % saving polar coords
    theta = polar_coords(electrode_index,2); 
    rho   = polar_coords(electrode_index,3); 
    
    % converting to cartesian coords
    [x,y] = pol2cart(theta,rho); 
    
    % appending row to our array 
    cart_coords = [cart_coords; electrode_index,x,y,]; 
end


%% Using Page 280 of Analyzing Neural Time Series Data 
numElectrodes = numel(EEG.chanlocs.X);


% The final laplacian voltage
ViLaplacian = zeros(numElectrodes); 

for eindx=1:numElectrodes

    % need to find the surrounding electrode within the boundaries 
    x = cart_coords(eindx,2); 
    y = cart_coords(eindx,3);

    electrodes_within_radius = []; 
    % Storing electrode indexes that are within our given radius
    for radIndx=1:length(cart_coords)
        if (abs(x - cart_coords(radIndx,2)) > radius(1) && abs(x-cart_coords(radIndx,3)) < radius(2))
            electrodes_within_radius = [electrodes_within_radius,radIndx]; 
        end
    end


    % taking all data at one electrode 
    Vi = EEG.data(x,:,:);
end


% Extract epochs from 0 to 800 ms. Compute the event-related potentials (ERPs) at each
% electrode. Get the average over all the trials. Select nine time points at which to show
% topographical plots (e.g., 0 to 800 ms in 100ms steps). In one figure, make a series of
% topographical plots at these time points. To increase the signal-to-noise ratio, make each plot
% show the average of activity from 20 ms before until 20 ms after each time point. For example,
% the topographical plot from 200 ms should show average activity from 180 ms until 220 ms.
% Indicate the center time point in a title on each subplot.

epochWindow = [-22 900]; 

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
num_rows = 2;
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

    % if you dont want the averages -- compute 
    % erp(:,timepointindex)
   
    % Create a subplot
    subplot(num_rows, num_cols, i);
    
    % Plot the topographic map
    topoplot(double(averages),elocFile,'eeg'); 
    title([num2str(time_point) 'ms']);
    
    c= colorbar; 
    colorbar;
    caxis([-10 10]);
end

colormap('jet')



