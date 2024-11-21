% eeg_topoplot displays a topoplot using 
% an instance of EEG
function eeg_topoplot(EEG, time_point_in_ms)

    % Assuming EEG contains your EEG data
    data_to_plot = EEG.data; % Replace with your actual EEG data

    % displaying the topoplot
    figure;
    topoplot(data_to_plot(:, time_point_in_ms), EEG.chanlocs); % EEG.chanlocs contains electrode positions
    title(['Topography at Time Point ', num2str(time_point_in_ms)]);
    colorbar; % Add a colorbar to the plot

end