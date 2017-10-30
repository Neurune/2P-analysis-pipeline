% Load in raw fluorescence signal extracted using MIJ

BaseDirectory = uigetdir('');
cd(BaseDirectory);
[filename] = uigetfile('*.mat','Select the MATLAB code file’);
load(filename); Grab the cells.mat file containing raw fluorescence for all ROIs

%% Compute moving baseline

cells.rawF = cells.rawF';

for i=1:size(cells.rawF,1);
    cells.baseline(i,:) = prctileFilt(cells.rawF(i,:),800); %10th percentile filter 
    disp({'Computing baseline for cell', num2str(i)});
end

% Sanitify check
hFig = figure;
for i=1:size(cells.rawF,1);
    hold on; plot(cells.rawF(i,:),'k')
    plot(cells.baseline(i,:),'r');
    %axis([-inf inf min(cells.rawF(:,i)) max(cells.rawF(:,i))*1.3 ])
    title(['ROI number ', num2str(i)]);
    ylabel('Raw fluorescence'); 
    xlabel('Time (in frames)');
    goodplot;
    pause;
    clf;
end
close(hFig);

%% Find dF/F using the moving baseline
cells.df = (cells.rawF-cells.baseline)./cells.baseline; 

time_s = linspace(0,length(cells.df)/30,length(cells.df));
% Sanitify check
hFig = figure;
for i=1:size(cells.df,1);
    hold on; plot(time_s,cells.df(i,:),'k');
    %hold on; plot(time_s,smooth(cells.df(i,:),30),'r');
    axis([-inf inf min(cells.df(i,:)) 0.8 ])
    title(['ROI number ', num2str(i)]);
    ylabel('dF/F0'); 
    xlabel('Time [seconds]');
    goodplot;
    pause;
    clf;
end
close(hFig);

%% Save cells variables
uisave('cells','cells'); %Will prompt you to confirm where to save the cell variable

clc; clear all; close all;