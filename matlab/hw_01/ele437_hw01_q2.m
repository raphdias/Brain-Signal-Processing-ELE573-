% Question B.) 

% Loop through each electrode and find the peak time of the ERP between 
% 100 and 800 ms. Store these peak times in a separate variable and then 
% make a topographical plot of the peak times (that is, the topographical 
% map will illustrate times in milliseconds, not activity at peak times). 
% Include a color bar in the figure and make sure to show times in 
% milliseconds from time 0 (not, for example, time indices instead of 
% milliseconds). What areas of the scalp show the earliest and the latest 
% peak responses to the stimulus within this window?

%% Start

% --- After running ele437_hw01_q1.m ---
% You will have a variable stored called 'erp'
% This variable we contain a 3D matrix with each electrodes time-averaged 
% value. 

% save the size of the erp
erp_data_size = size(erp); 

% save the time indexes to loop through
start_time = 100; % ms
end_time   = 800; % ms

start_index = round(((start_time - epochedTimeMatrix(1)) / 1000) * EEG.srate) + 1; 
end_index   = round(((end_time - epochedTimeMatrix(1))   / 1000) * EEG.srate) + 1; 

% creating a new time array to store times between 100ms and 800ms 
timeMatrixFor100to800ms = epochedTimeMatrix(start_index:end_index); 

% loop through each electrode row and store the maximum value 
% and peak times 
% between the start_index and the end_index
max_values = []; 
peak_times = []; 

for electrode_index = 1:erp_data_size(1)

    % isolating the data we want analyze for an electrode
    electrode_data = erp(electrode_index,start_index:end_index); 

    % find the max value of electrode_data 
    [max_val, max_index]  = max(electrode_data(:));

    % append max_value and max_time to a list 
    max_values(length(max_values) + 1,1) = max_val;
    peak_times(length(peak_times) + 1, 1) = timeMatrixFor100to800ms(max_index); 
end

% Change the color bar scale 


% topoplot of peak times
topoplot(peak_times,elocFile,'EEG', 'electrodes','labels'); 
title('Peak times of erp between 100 and 800 ms');
c = colorbar;
caxis([0 600]);
caxis([0 600]);
c.Label.String = 'Time in ms';
colormap('Jet')





