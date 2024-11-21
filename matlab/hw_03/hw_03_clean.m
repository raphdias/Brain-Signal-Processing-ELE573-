% Homework 3

close all; 
clear; 

%% Custom Params

% plotting 
plotMorlets = true;

% setting working directory 
cd('/Users/raph/Desktop/Fall 2024/ele537/matlab/hw_03')
addpath('data'); 

load('sampleEEGdata.mat'); 
elocFile = 'eloc64C2.txt'; 


%% Question A 

morletFreqz = linspace(2, 30, 5); % cmw frequencies 
morN = length(morletFreqz); % number of freqz
cycles = 4; 
freqScale = 1;
Fs = EEG.srate;
dt = 1/Fs; 

cmwFamily = cell(1,length(morletFreqz));  % empty cell to store cmw

% creating and storing cmw
for i = 1:morN
    
    f = morletFreqz(i); 
    t = -1:dt:1;
    sigma = cycles / (2 * pi * f);

    sine = 1/sqrt(sigma*pi^0.5)*exp(1i * 2 * pi * f * t); 
    gauss = exp(-t.^2 / (2 * sigma^2)); 

    cmwFamily{i} = freqScale.*sine.*gauss; 

end

if plotMorlets
    figure(); 
    for i = 1:morN
        subplot(morN,1,i)
        plot(t,cmwFamily{i})
        set(gca,'ytick',[])
        title([num2str(morletFreqz(i)) ' Hz'])
    end
    sgtitle('Complex Morlet Wavelets')
end

clear i f sigma freqScale sine gauss 

%% Question B

eeg_data = EEG.data(:,1:length(t),1); % only taking length of cmw
[num_electrode, num_samples] = size(eeg_data); 

convolved_data = cell(1, morN); % Create an empty cell array

% computing the convolution in freq domain
for i = 1:morN
    cmw_fft = fft(cmwFamily{i});
    temp_conv = []; 
    for k = 1:num_electrode
        temp_conv = [temp_conv;ifft(cmw_fft.*fft(eeg_data(k,:)))]; % convolving eeg_fft with each cmw
    end
    convolved_data{i} = temp_conv;
end

clear i k cmw_fft temp_conv

%% Question C 

powerPhaseMatrix = zeros(width(t), 5, 64, 2);  % create 4D-Matrix 

for i = 1:morN
    for k = 1:num_electrode
        % take power of the electrode 
        c_data = convolved_data{i}(k,:); 
        c_power = abs(c_data).^2; 
        c_angle = angle(c_data); 

        powerPhaseMatrix(:,i,k,1) = c_power; 
        powerPhaseMatrix(:,i,k,2) = c_angle; 
    end
end

%% Part D and E - Repeat step (d) for activity at 360 ms, and 650 ms.

power_and_phase_plots(180, powerPhaseMatrix, EEG, elocFile,morletFreqz);
power_and_phase_plots(360, powerPhaseMatrix, EEG, elocFile,morletFreqz);
power_and_phase_plots(650, powerPhaseMatrix, EEG, elocFile,morletFreqz);

%%

function power_and_phase_plots(time_ms, powerPhaseMatrix, EEG, elocFile,morletFreqz)


    indexms = round((time_ms/1000 - (EEG.times(1))/1000) * (EEG.srate)) + 1;
    
    columns = 5;  % frequncies 
    rows = 2;     % power/phase
    
    % calculating global mins and maxes for power
    power_plot = double(squeeze(powerPhaseMatrix(indexms,:,:,1))); 
    
    % global min and max of power
    power_log_max = log(max(power_plot(:)));
    power_log_min = log(min(power_plot(:)));
    
    % Create subplots and store their handles in a cell array
    subplotHandles = cell(rows, columns);
    
    figure(); 
    % plotting power
    for i = 1 : columns
    
        ax(1,i) = subplot(rows,columns,i); 
        power = squeeze(powerPhaseMatrix(indexms, i, :, 1)); 
        topoplot(log(power),elocFile, 'EEG'); 
        caxis([power_log_min power_log_max])
        title(['Frequency - ' num2str(morletFreqz(i))]);
        ax(2,i+columns) = subplot(rows,columns,i+columns); 
        topoplot(squeeze(powerPhaseMatrix(indexms, i, :, 2)),elocFile, 'EEG');
    
    end
    
    nColors = 256; 
    exp_colormap = jet(nColors); 
    colormap(exp_colormap);
    
    
    c1 = colormap(ax(1,columns),'jet');
    pos = get(subplot(rows,columns,columns),'Position');
    h = colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)  pos(3)/10  pos(4)]);
    
    c2 = colormap('jet');
    pos = get(subplot(rows,columns,columns*2),'Position'); 
    h = colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)  pos(3)/10  pos(4)]);

    sgtitle(['Colormaps at ' num2str(time_ms) ' ms'])
end




 
