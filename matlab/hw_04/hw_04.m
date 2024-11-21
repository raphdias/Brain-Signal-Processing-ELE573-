% HW 04
clear; close all; clc; 

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
time_window_after = [100 602];

% Extract data for the two time windows
indices_before = find(EEG.times >= time_window_before(1) & EEG.times <= time_window_before(2));
eeg_data_before = EEG.data(:,indices_before,:);

indices_after = find(EEG.times >= time_window_after(1) & EEG.times <= time_window_after(2));
eeg_data_after = EEG.data(:,indices_after,:);


% Center data for the time window
eeg_pre = eeg_data_before - mean(eeg_data_before, 2);
eeg_post = eeg_data_after - mean(eeg_data_after, 2);

clear eeg_data_after eeg_data_before indices_after ; 

%% Covariance

covar= zeros(2,EEG.nbchan, EEG.nbchan, EEG.trials);
% for each trial, obtain sample cov
for triali = 1:EEG.trials
    eeg_b = squeeze(eeg_pre(:,:,triali)); 
    eeg_b = (eeg_b*eeg_b')./(length(eeg_b)-1);
    covar(1,:,:,triali) = eeg_b;

    eeg_f = squeeze(eeg_post(:,:,triali)); 
    eeg_f = (eeg_f*eeg_f')./(length(eeg_f)-1);
    covar(2,:,:,triali) = eeg_f;
end

clear eeg_b eeg_f triali

cov_erp = mean(covar,4);

subplot(1,2,1); 
imagesc(squeeze(cov_erp(1,:,:)))
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[20 150])
title('Average Covariance of EEG over trials - Before')

subplot(1,2,2); 
imagesc(squeeze(cov_erp(2,:,:)))
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[20 150])
title('Average Covariance of EEG over trials - After')


%% Eigen 

covar_pre = squeeze(cov_erp(1,:,:)); 
covar_post = squeeze(cov_erp(2,:,:)); 

% principle components analysis via eigenvalue decomposition
[pc_pre,eigvals_pre] = eig(covar_pre);
[pc_post,eigvals_post] = eig(covar_post);


% components are listed in increasing order, and converted here to descending order for convenience
pc_pre      = pc_pre(:,end:-1:1);
eigvals_pre = diag(eigvals_pre);
eigvals_pre = 100*eigvals_pre(end:-1:1)./sum(eigvals_pre); % convert to percent change
% eigvals_pre = eigvals_pre(end:-1:1); % convert to percent change

pc_post      = pc_post(:,end:-1:1);
eigvals_post = diag(eigvals_post);
eigvals_post = 100*eigvals_post(end:-1:1)./sum(eigvals_post); % convert to percent change
% eigvals_post = eigvals_post(end:-1:1); % convert to percent change


%% plot topomap and time course for first 4 components

figuretxt = 'Four Principal Components of pre and post stimulus'; 
F = figure('Name',figuretxt,'NumberTitle','off');
T=tiledlayout(F,'Flow','TileSpacing','compact','Padding','none');
for i=1:4

    t = tiledlayout(T,'flow','TileSpacing','tight','Padding','none');
    t.Layout.Tile = i;

    nexttile(t);
    topoplot(double(abs(pc_pre(:,i))),elec_file_path,'eeg');
    caxis([0 0.2]);
    title(['Window= [-500 0] ms ', 'eigval=' num2str(eigvals_pre(i)) '%'])

    nexttile(t);
    topoplot(double(abs(pc_post(:,i))),elec_file_path,'eeg');
    caxis([0 0.2]);
    title(['Window= [100 600] ms ', 'eigval=' num2str(eigvals_post(i)) '%'])

    
    % PC time course for each trial, then average together
    pctimes = zeros(1,EEG.pnts);
    for triali=1:EEG.trials
        eeg = bsxfun(@minus,squeeze(EEG.data(:,:,triali)),squeeze(mean(EEG.data(:,:,triali),2)));
        pctimes = pctimes + pc_pre(:,i)'*eeg;
    end

    nexttile(t);
    plot(EEG.times,pctimes./EEG.trials)
    ylabel('voltage')
    set(gca, 'ylim', [-40 40])
    set(gca, 'xlim', [-1000 1500])
    

    % PC time course for each trial, then average together
    pctimes = zeros(1,EEG.pnts);
    for triali=1:EEG.trials
        eeg = bsxfun(@minus,squeeze(EEG.data(:,:,triali)),squeeze(mean(EEG.data(:,:,triali),2)));
        pctimes = pctimes + pc_post(:,i)'*eeg;
    end

    nexttile(t);
    plot(EEG.times,pctimes./EEG.trials)
    set(gca, 'ylim', [-40 40])
    set(gca, 'xlim', [-1000 1500])
    sgtitle(['PC # ' num2str(i)])
end

axes(T,'visible','off');
cb = colorbar;
caxis([0 0.2]);
colormap('Jet')
cb.Layout.Tile = 'east';
cb.Label.FontSize = 25;
cb.Label.String = 'Eigen Vector';
xlabel(T,'Time (ms)','fontsize',25)
title(T,figuretxt,'fontsize',25)
% set(F,'WindowState','normal','WindowStyle','docked')
