close all
clear 
clc

%% Load data
Fs = 1000; %[Hz]  

dynamometer = load('dynamometer.mat');

[num_elements,num_channels] = size(dynamometer.data);
t = (0:num_elements-1)/Fs;
figure('Name','Dynamometer')
for i = 1:num_channels
    subplot(num_channels,1,i)
    plot(t,dynamometer.data(:,i))
    title(dynamometer.labels(i,:))
    ylabel(dynamometer.units(i,:))
    xlabel('t [s]')
end

EMG = load('Electromyography.mat');


[num_elements,num_channels]=size(EMG.data);

t = [0:num_elements-1]/Fs;
figure('Name','EMG')
for i = 1:num_channels
    subplot(num_channels,1,i)
    plot(t,EMG.data(:,i))
    title(EMG.labels(i,:))
    ylabel(EMG.units(i,:))
    xlabel('t [s]')
end


%%
sig = EMG.data(:,1);    % Work on a single channel of EMG

figure
plot(t,sig)
title('Observed signal')

%% Band-pass filtering

BP = load('BP_EMG1','-ascii');
s=filter(BP,1,sig);

figure
N=length(sig);
X=fftshift(abs(fft(s)))/N;
f=((0:N-1)-floor(N/2))*Fs/N;
plot(f,X)
title('Frequency spectrum of the filtered signal')
figure
subplot(2,1,1)
plot(t,s)
title('Filtered signal')
subplot(2,1,2) 
x=abs(s);
plot(t,x)
title('Rectified signal')

%% ENVELOPES

N = 0.5*Fs; % Half-second window  

%% IEMG

b=zeros(1,N+1);
b(1)=1;
b(N+1)=-1;
a=[1,-1];

x_IEMG=filter(b,a,x);
figure
subplot(2,1,1)
plot(t,x)
hold on
plot(t,x_IEMG)
legend('Rectified signal','IEMG','Location','best')
subplot(2,1,2)
plot(t,x_IEMG)
title('IEMG')

%% ARV

b=zeros(1,N+1);
b(1)=1;
b(N+1)=-1;
a=N*[1,-1];
x_ARV=filter(b,a,x);


figure
subplot(2,1,1)
plot(t,x)
hold on
plot(t,x_ARV)
legend('Rectified signal','ARV','Location','best')
subplot(2,1,2)
plot(t,x_ARV)
title('ARV')

%% RMS

x_RMS = s.^2;
b=zeros(1,N+1);
b(1)=1;
b(N+1)=-1;
a=N*[1,-1];
x_RMS=filter(b,a,x_RMS);

figure
subplot(2,1,1)
plot(t,x)
hold on
plot(t,x_RMS)
legend('Rectified signal','RMS','Location','best')
subplot(2,1,2)
plot(t,x_RMS)
title('RMS')

%% Crossing rate using threshold previously computed
p=length(s);
vector=zeros(p+N,1);
vector(1:p) = s;
CR=zeros(p,1);
rms =  s(5*Fs:6.5*Fs).^2;
rms = sqrt(mean(rms));

for i=0:(length(vector)-N-1)
    
    CR(i+1,1) = crossing_rate(vector((i+1):i+N),2*rms);
    
end
figure
% time = length(CR)/fc;
subplot(2,1,1)
plot(t,s)
hold on
plot(t,CR)
legend('Rectified signal','Crossing rate','Location','best')

subplot(2,1,2)
plot(t,CR)
title('Crossing rate')


%% Adapt amplitude range of the envelope 0 - 100
figure

x_RMS = x_RMS - min(x_RMS);
x_RMS = x_RMS.*100./max(x_RMS);
plot(t,x_RMS)
hold on

x_IEMG = x_IEMG - min(x_IEMG);
x_IEMG = x_IEMG.*100./max(x_IEMG);
plot(t,x_IEMG)
hold on

x_ARV = x_ARV - min(x_ARV);
x_ARV = x_ARV.*100./max(x_ARV);
plot(t,x_ARV)
hold on

CR = CR - min(CR);
CR = CR.*100./max(CR);
plot(t,CR)

legend('RMS', 'IEMG', 'ARV', 'Crossing rate','Location','best')

title('Normalized envelopes')