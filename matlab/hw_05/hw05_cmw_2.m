%% Matlab Code
% BME 473/ELE573 Course Assignment3
clear all
load sampleEEGdata.mat
data=EEG.data;
srate = EEG.srate; % sampling rate, sometimes abbreviated Fs \
% the time in EEG is in ms in the formula we must use it in seconds to have
the correct wave
time = EEG.times/1000;
% a) Create a family of Morlet wavelets from 2 Hz to 30 Hz in five steps
(fixed cycle=4)
%we break this into 3 steps all of them are going to be implemented in a
%”for” loop
% Step 1) Create the Complex sin wave corresponding to each frequency f
% Step 2) Create the Gaussian wave
% Step 3) Multiply the Complex sine Wave and the Gaussian window (using
% element wise multiplication)
%you will need to define the time vector:
t = linspace(-1,1,2*EEG.srate);
% or:
% t=-1:1/EEG.srate:1-(1/EEG.srate);% we use this -1/EEG.srate trick to get
exactly 512 samples for 2 seconds
% you will need to define your frequemncy vector:
% you can set it manually freq_vector=[2 9 16 23 30];
% or you can define it as follows
max_freq=30;
min_freq=2;
step=(max_freq-min_freq)/4; % 30 Hz is included
freq_vector=min_freq:step:max_freq;
Number_of_freq=length(freq_vector);
Number_of_ch=64;
n_cycles=4;
% you can also use linspace as follows:
% freq_vector=round(linspace(min_freq, max_freq, 5));
% Then you will loop on each frequency in order to create your complex sine
% wave and Gaussian and multiply them together and store the result in a
% matrix "Morlet"
for i=1:Number_of_freq

 f=freq_vector(i);
 s=n_cycles/(2*pi*f);
 A=1/sqrt(s*pi^0.5);%complex sine wave amplitude
 ComplexSineWave=A*exp(1i*2*pi*f*t);
 %Gaussian wave , 4 cycles
 GW=exp(-1*t.^2/(2*s.^2));
 Morlet(i,:)=ComplexSineWave.*GW;%creating our wavelet by element wise
multiplication
end
%if you want to plot in order to check your wavelets, plot the real part
% for i=1:5
% plot(t,real(Morlet(i,:)));
% hold on
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b) Convolve each wavelet with EEG data from all electrodes and from only
first trial
% hint: implement the convolution using the multiplication in the frequency
% domain (more efficient in Matlab than time domain convolution)
clear TrialData SizeData SizeMorlet SizeConv SizeHalfMorlet FFTData FFTMorlet
clear Prod CONV CONVData
% extract the data of a single trial of all the channels
TrialData=data(:,:,1);
% The Size = the number of samples
SizeData = size(TrialData,2);
SizeMorlet = size(Morlet,2);
% prepare the size of the convolution result
SizeConv = SizeData+SizeMorlet-1;
% After the convolution in the frequency domain, we need the final size to
% be equal to the size of the data so we will cut half the size of the
% Morlet wavelet from the begining and end of the data
SizeHalfMorlet = floor(SizeMorlet/2);
% loop on the channels and the frequencies
for ch=1:Number_of_ch
 for i=1:Number_of_freq
 %then FFT
 FFTData = fft(TrialData(ch,:),SizeConv);%the size of Conv.
 FFTMorlet = fft(Morlet(i,:),SizeConv);
 % Normalize the wavelet in the frequency domain
 FFTMorlet = FFTMorlet./max(FFTMorlet);
 % Multiply then inverse FFT

 Prod=FFTData.*FFTMorlet;
 CONV = ifft(Prod);

 % cut from beginning and end
 CONV = CONV(SizeHalfMorlet:end-SizeHalfMorlet);
 %Store The Result of the complex wavelet convolution
 %if needed check https://www.mathworks.com/help/matlab/ref/colon.html
for the use of the column operator in indexing matrices
 CONVData(:,i,ch)=CONV;
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c) Extract power and phase from the result of the complex wavelet
%convolution and store in a time*frequency*channels*power/phase
clear M CWPower CWTheta CWPowerPhase plottimevec freqvec tplot PowerToPlot
clear PhaseToPlot vec
%Note that we did not use loops here as the abs and angle functions operate on
the
%whole matrix and returns a similar size matrix
Power=abs(CONVData).^2;
Phase=angle(CONVData);
% or:
%M=sqrt(real(CONVData).^2+imag(CONVData).^2);
%Power=M.^2;
% Phase=atan(imag(CONVData)./real(CONVData)); 
%Store in a a time*frequency*channels*power/phase Matrix as required
%Note the column operator transfer all the dimensions from Power and Phase
%as they are then we add the fourth dimension as we want to store the Power
%value in 1 and the Phase in 2
CWPowerPhase(:,:,:,1)=Power;
CWPowerPhase(:,:,:,2)=Phase;
%Using the column operator is equivalent to using the for loops on each
dimension as follows:
%the loop implementation in matlab is not efficient in terms of time. A
%very useful practice in Matlab programming is to use the functions and
%operators to operate directly on matrices instead of looping on dimensions
% for s=1:SizeData
% for i=1:Number_of_freq
% for ch=1:Number_of_ch
%
%
% CWPowerPhase(:,:,:,1)=Power;
% CWPowerPhase(:,:,:,2)=Phase;
%
% end
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%d) Make topographical plots of power and phase at 180 ms, 360ms, and 650ms
at all frequencies
%set the time vector, note that we add 1000 ms to each time as the
%EEG.times starts 1 sec before the stimulus
plottimevec=[1180 1360 1650]; timevec=[180 360 650];
%creating the 3 required plots in a for loops, we use the indexes of the
%loops in creating the figures
close all
for plotindex=1:3
 figure(plotindex)
 %calculate the samples corresponding for the time for plotting
 tplot=round(plottimevec(plotindex)/1000*srate);

 PowerToPlot=squeeze(CWPowerPhase(tplot,:,:,1));
 PhaseToPlot=squeeze(CWPowerPhase(tplot,:,:,2));
 for j=1:Number_of_freq
 subplot(2,5,j);
 vec=PowerToPlot(j,:);
% l=min(vec);
% u=max(vec);
 topoplot(double(vec),'eloc64C2.txt','eeg');
 titlestr=sprintf('Power at %d ms at %d
Hz',timevec(plotindex),freq_vector(j));
 title(titlestr, 'fontsize', 15); 
% caxis([l u]);
 colormap(jet);
 colorbar;
 end
 for j=1:5
 subplot(2,5,j+5);
 vec=PhaseToPlot(j,:);
% l=min(vec);
% u=max(vec);

 topoplot(double(vec),'eloc64C2.txt','eeg');
 titlestr=sprintf('Phase at %d ms at %d Hz',timevec(plotindex),
freq_vector(j));
 title(titlestr, 'fontsize', 15);
% caxis([l u]);
 colormap(jet);
 colorbar;
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part g)
%Consider Frequency Range from 2 to 40 Hz varying in 2 steps
%and varying cycles [3-10] %code is repeated for fixed cycle too
clear step index Morlet s A GW ComplexSineWave
step=2;
index=1;
% define the vector for the variable number of cycles
cyclevec=[3 4 5 6 7 8 9 10];
%loop over frequencies from 2 to 40 hz with step of 2
 for f=2:step:40

 %define the wavelet using the fixed cycle uses n_cycles=4 like the previous
code
 Fs=4/(2*pi*f);
 FA=1/sqrt(Fs*pi^0.5);
 FComplexSineWave=FA*exp(1i*2*pi*f*t);
 %Gaussian wave , 4 cycles
 FGW=exp(-1*t.^2/(2*Fs.^2));
 FMorlet(index,:)=FComplexSineWave.*FGW;%creating our wavelet by
multiplication

 %define the wavelet using variable number of cycles
 %choose the number of cycles
 cindex=round(f/6)+1;
 s=cyclevec(cindex)/(2*pi*f);
 A=1/sqrt(s*pi^0.5);
 ComplexSineWave=A*exp(1i*2*pi*f*t);
 GW=exp(-1*t.^2/(2*s.^2));
 Morlet(index,:)=ComplexSineWave.*GW;%creating our wavelet by
multiplication
 index=index+1;
 end
clear OTrialData OFFTData OFFTMorlet Prod OCONV OCONVData
%Channel FCz only
OTrialData=squeeze(data(47,:,:));
OTrialData=OTrialData';%trials*samples
% first the sizes
SizeData = size(OTrialData,2);
SizeMorlet = size(Morlet,2);
SizeConv = SizeData+SizeMorlet-1; SizeHalfMorlet = floor(SizeMorlet/2);
for trial=1:99
 OFFTData = fft(OTrialData(trial,:),SizeConv);%the size of Conv
 for i=1:20 %number of frequencies , now it is 20
 %then FFT
 OFFTMorlet = fft(Morlet(i,:),SizeConv);
 %Normalization of the wavelet in the frequency domain
 OFFTMorlet= OFFTMorlet./max(OFFTMorlet);

 FOFFTMorlet = fft(FMorlet(i,:),SizeConv);%convolution with fixed
morlet
 %Normalization of the wavelet in the frequency domain
 FOFFTMorlet=FOFFTMorlet./max(FOFFTMorlet);
 % Multiply then inverse
 Prod=OFFTData.*OFFTMorlet;
 OCONV = ifft(Prod);
 FProd=OFFTData.*FOFFTMorlet;
 FOCONV = ifft(FProd);

 % cut from beginning and end
 OCONV = OCONV(SizeHalfMorlet:end-SizeHalfMorlet);
 OCONVData(trial,i,:)=OCONV;%contains for each of the 99 trials of FCz
the convolution for each of the 20 frequencies

 FOCONV = FOCONV(SizeHalfMorlet:end-SizeHalfMorlet);
 FOCONVData(trial,i,:)=FOCONV;
 end
end

%get the power for each trial
OCWPower=abs(OCONVData).^2;
FOCWPower=abs(FOCONVData).^2;
% calculate the mean of the power over all the trials
OMeanCWPower=squeeze(mean(OCWPower,1));
FOMeanCWPower=squeeze(mean(FOCWPower,1));

%Apply Baseline Correction
Baseline=[-500 -200];
%Claculate the samples corresponding to the Baseline period
BaselineInd=[round((Baseline(1)/1000)*EEG.srate)+EEG.srate
round((Baseline(2)/1000)*EEG.srate)+EEG.srate];
BaselinePower=mean(OMeanCWPower(:,BaselineInd(1): BaselineInd(2)),2);
FBaselinePower=mean(FOMeanCWPower(:,BaselineInd(1): BaselineInd(2)),2);
%Now every frequency is divided by the avg corresponding to this frequency
for i=1:20
 OMeanCWPower(i,:)=OMeanCWPower(i,:)./BaselinePower(i);
 FOMeanCWPower(i,:)=FOMeanCWPower(i,:)./FBaselinePower(i);
end
%%Plot Time Frequency MAP 
%x: time [-200 1000 ms], y: frequency [2 40hz] and color 10log10(power)
OMeanCWPower=OMeanCWPower';
FOMeanCWPower=FOMeanCWPower';
% Plotting
figure(7),
subplot(1,2,1)
xaxis=EEG.times(206:512);
yaxis=linspace(2,40,20);
MapToPlot=OMeanCWPower(206:513,:);
MapToPlot=MapToPlot';
imagesc(xaxis,yaxis,10*log10(MapToPlot));
caxis([-3 3]);
colorbar
colormap('jet');
set(gca,'YDir','normal')
title('Time-Frequence Map: variable cycles')
xlabel('Time ms')
ylabel('Frequency Hz')

subplot(1,2,2)
FMapToPlot=FOMeanCWPower(206:512,:);
FMapToPlot=FMapToPlot';
imagesc(xaxis,yaxis,10*log10(FMapToPlot));
caxis([-3 3]);
colorbar
colormap('jet');
set(gca,'YDir','normal')
title('Time-Frequence Map: fixed (4) cycles')
xlabel('Time ms')
ylabel('Frequency Hz')