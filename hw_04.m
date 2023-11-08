% HW 04

%% Imports 

data_folder_path = '/Users/raph/Desktop/ele537/matlab/data_sets'; 

load(fullfile(data_folder_path,'sampleEEGdata.mat'));
elec_file_path = fullfile(data_folder_path,'eloc64C2.txt');
 
% Open the file for reading
fileID = fopen(elec_file_path, 'r');
 
% Read the data using textscan
data_elec_file = textscan(fileID, '%d %f %f %s', 'Delimiter', '\t');

clear fileID

%% PCA broadband

Fs = 256; % Hz

% Define the time windows in milliseconds
time_window_before = [-500 0];
time_window_after = [100 600];

% Extract data for the two time windows
indices_before = find(EEG.times >= time_window_before(1) & EEG.times <= time_window_before(2));
eeg_data_before = EEG.data(:,indices_before,:);

indices_after = find(EEG.times >= time_window_after(1) & EEG.times <= time_window_after(2));
eeg_data_after = EEG.data(:,indices_after,:);


% Center data for the time window
eeg_data_before_centered = eeg_data_before - mean(eeg_data_before, 2);
eeg_data_after_centered = eeg_data_after - mean(eeg_data_after, 2);

%% Covariance


cov_matrices_before = zeros(64, 64, width(eeg_data_before_centered));
cov_matrices_after = zeros(64, 64, width(eeg_data_after_centered));


% Calculate the sample covariance matrices for the 'before' time window
for i = 1:EEG.trials
    eeg_data_before = eeg_data_before_centered(:,:,i);
    eeg_data_after  = eeg_data_after_centered(:,:,i);

    cov_matrices_before(:,:,i) = (eeg_data_before*eeg_data_before')/(width(eeg_data_before)-1); 
    cov_matrices_after(:,:,i) = (eeg_data_after*eeg_data_after')/(width(eeg_data_after)-1); 
end

clear i

% average of covariances matrix
cov_before_avg = mean(cov_matrices_before,3);
cov_after_avg = mean(cov_matrices_after,3);


% Eigenvalue decomposition for the 'before' time window
[eigenvectors_before, eigenvalues_before] = eig(cov_before_avg);
[eigenvalues_before, order_before] = sort(diag(eigenvalues_before), 'descend');
eigenvectors_before = eigenvectors_before(:, order_before);

% Eigenvalue decomposition for the 'after' time window
[eigenvectors_after, eigenvalues_after] = eig(cov_after_avg);
[eigenvalues_after, order_after] = sort(diag(eigenvalues_after), 'descend');
eigenvectors_after = eigenvectors_after(:, order_after);