% eeg_describe gives you information on what your eeg data
% contains 
function eeg_describe(EEG)

    % Check if the EEG structure is valid
    if isempty(EEG)
        error('EEG structure is empty. Please provide a valid EEG instance.');
    end

    % Display EEG properties
    fprintf('*** EEG Information ***\n');
    fprintf('Number of Channels: %d\n', EEG.nbchan);
    fprintf('Number of Time Points: %d\n', EEG.pnts);
    fprintf('Number of Trials: %d\n', EEG.trials);

    % Display additional information if available
    if isfield(EEG, 'srate')
        fprintf('Sampling Rate: %d Hz\n', EEG.srate);
    end

    if isfield(EEG, 'subject')
        fprintf('Subject: %s\n', EEG.subject);
    end

end