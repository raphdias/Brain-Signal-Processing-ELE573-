newEEGDataStorage = zeros(size(EEG.data)); 

for electrode = 1:EEG.nbchan

    x = cart_coords(electrode,2); 
    y = cart_coords(electrode,3); 
    
    e_within_rad = []; 
    for radIndx=1:numel(cart_coords(:,1))
    
       % compute the distance between the two electrodes and see whether they
       % are in the defined range 
       d = sqrt ((cart_coords(radIndx,2) - x)^2 + (cart_coords(radIndx,3) - y)^2); 
    
       % if our distance is within the radius, we can mark the electrode
       if (d > radius(1) && d < radius(2))
           e_within_rad = [e_within_rad; radIndx d];
       end
    end
    
    % once we have collected all of our electrodes within the radius
    % we can compute g
    gNetElectrode = []; 
    
    % THIS IS WEIGHTING FORMULA IS PROBABLY WRONG 
    numElectrodesWithinRad = numel(e_within_rad(:,1)); 

    distanceSum = 0; 
    
    for sumIndx=1:numElectrodesWithinRad
        distanceSum = distanceSum + 1/e_within_rad(sumIndx,2); 
    end

    for gIndex=1:numElectrodesWithinRad

        gijSingle = (1/e_within_rad(gIndex,2))/(distanceSum);

        gNetElectrode = [gNetElectrode; gIndex gijSingle]; 
    end
    
    
    % Calculating the Voltage depending on the finite difference method
    electrode_of_interest = EEG.data(electrode,:,:); 
    
    % loop through each voltage point, at each time
    % theoretical loop thorugh each electrode
    
    [e, volts, trials] = size(electrode_of_interest);
    for timeIndex = 1:trials 
        for posIndex = 1:volts
            voltOnPrimaryElectrode=electrode_of_interest(1,posIndex,timeIndex);
    
            % caluclating net gij * Vj
            gijSum = 0; 
            for surroundElectrodes=1:numElectrodesWithinRad
                gijSum = gijSum + gNetElectrode(surroundElectrodes,2) * EEG.data(gNetElectrode(surroundElectrodes,1),posIndex,timeIndex);
            end
            newEEGDataStorage(electrode, posIndex, timeIndex) = voltOnPrimaryElectrode - gijSum; 
        end
    end
end

EEG.data = newEEGDataStorage;



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
sgtitle('Large Laplacian Applied To Data');
colormap('jet')
