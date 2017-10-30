%% Script for getting raw fluorescence signal out of a set of ROIs already defined and saved. 

%  The code used MIJ to interact with FIJI/ImageJ via MATLAB. A set of ROIs
%  (already saved in a .zip file is loaded in at the raw fluorescence for
%  all ROIs is extracted and saved in a cell.

%  The code is tested on a Mac in MATLAB 2015a 
%  It's a good idea to increase the JAVA Heaps size in MATLAB before loading large stacks
%
%  Updated on: October 30, 2017 (Rune Rasmussen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start on a fresh
clear all; clc; 

%% Change the Java Heap Memory if MIJ runs out of memory when loading in stacks
oldMaxHeapSize = com.mathworks.services.Prefs.getIntegerPref('JavaMemHeapMax');  
com.mathworks.services.Prefs.setIntegerPref('JavaMemHeapMax', 8000); % Set to 8000 MB
%quit(); %Need to re-start MATLAB to take action

%% Set the current directory (cd) to where the registered image stack and ROI set is located

%The code way
    %BaseDirectory = '/Users/RuneRasmussen/Dropbox/PhD - Yonehara Lab/Calcium imaging analysis/Documentation/Test Image Stack';
    %cd(BaseDirectory);

%The click way - find the folder containing the registered image stack and the ROI set
    BaseDirectory = uigetdir('');
    cd(BaseDirectory);

%Pre-establish cell variable
cells = []; %This is the variable the everything will be saved into

%% Opening FIJI/ImageJ in MATLAB (MIJI) 
addpath('/Applications/Fiji.app/scripts'); %Make sure to add the FIJI scripts (where MIJ.m file is located) folder to your path
javaaddpath '/Applications/MATLAB_R2015b/java/jar/mij.jar' %MacBook path. This is not always needed, but if MATLAB 'forget's where this file is MIJ cannot run
%javaaddpath '/Applications/MATLAB_R2016a.app/java/jar/ij.jar' % iMac path. This is not always needed, but if MATLAB 'forget's where this file is MIJ cannot run


Miji; %This will open FIJI/ImageJ via MATLAB
import ij.*; %Some code that helps - don't ask me :) 
import ij.IJ.*; %Some code that helps - don't ask me :) 

%% Open registered image stack through with MIJ
MIJ.run('Open...'); %Choose the registered image stack
MIJ.run('Out [-]'); MIJ.run('Out [-]'); %Will decrease the size of the stack for convinience 
filename =  char(MIJ.getCurrentTitle); filename = filename(1:end-4); %Get rid of '.tif' 

%% Make average intensity Z-projection image for nice looking image

%Average Intensity projection image
    %MIJ.run('Z Project...', 'projection=[Average Intensity]');
    %cells.avgProject = MIJ.getCurrentImage; %Saves into the structure 'cells'
    %pause(3) %Allows you time to see the Z-projection image
    %MIJ.run('Close','');  
%Standard deviation intensity projection image
    %MIJ.run('Z Project...', 'projection=[Standard Deviation]');
    %cells.stdProject = MIJ.getCurrentImage; %Saves into the structure 'cells'
    %pause(3) %Allows you time to see the Z-projection image
    %MIJ.run('Close','');

%% Load in the ROI set into the ROI manager 
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
count = RC.getCount();
    if RC.getCount~=0 %Clears any ROIs already in memory
        RC.runCommand('Deselect');
        RC.runCommand('Delete');
    end
MIJ.run('Open...'); %Load the ROI set (i.e. the zip file
RC.runCommand('Show All'); %Will show all the ROI's on the stack and add a number to them

%% Load in the XY positions of all ROIs into MATLAB for later usage
import ij.*; import ij.IJ.*; RM = ij.plugin.frame.RoiManager(); RC = RM.getInstance();
cellNumber = RC.getCount(); %The number of cells (i.e. ROIs)

%Gets the position of each ROI, for making maps of ROI colour-coded for OS or DS later 
for i = 1:cellNumber
RC.select(i-1)
cells.loc(i).name = char(RC.getName(i-1));
tempString = char(cells.loc(i).name);
if length(tempString)>9
    [~,temp] = strtok(tempString,'-');
    [tempLoc1, tempLoc2] = strtok(temp,'-');
else
    [tempLoc1,tempLoc2] = strtok(tempString,'-');
end;
cells.xPos(i) =  str2num(tempLoc1)';
cells.yPos(i) =  -str2num(tempLoc2)';
end
clear tempLoc1; clear tempLoc2; clear tempString; 

%Sanity check if ROI position is correctly loaded in
% figure(1); c = linspace(1,10,cellNumber);
% scatter(cells.yPos, cells.xPos,120,c, 'filled'); axis([0 Xsize 0 Ysize]); 
% axis square; set(gca,'Ydir','reverse'); goodplot; box on; 
% pause(5); close

%% Get the time-series fluorescence intensity profiles (i.e. raw fluorescence) for all ROIs
RM = ij.plugin.frame.RoiManager(); 
RC = RM.getInstance();
    RC.runCommand('Deselect');
    RC.runCommand('Multi Measure'); %A dialog menu will pop up here, click 'OK' and a results table should be opened) !!!!!!

Headers = MIJ.getListColumns; %Gets the name of each column, corresponding to individual ROIs
%Now get the values for all ROI across N images
for i = 1:length(Headers)    
    cells.rawF(:,i) = MIJ.getColumn(Headers(i)); %Saves raw F values in cells structure
end
MIJ.run('Close',''); %Will close the result table
clear Headers; %Clears header variable

%Clear ROI manager and close it
count = RC.getCount();
if RC.getCount~=0 %Clears any ROIs already in memory
       RC.runCommand('Deselect');
       RC.runCommand('Delete');
end
MIJ.run('Close'); clear RC; clear RM;

%% Close all windows and close connection to MIJ (this sort of works)
MIJ.closeAllWindows; MIJ.exit;

%% Sanity checks

cellNumber = size(cells.rawF(:,:),2);
%Sanity check - Detailed: Toggle through every individual ROIs responses
figure(3);
for i=1:cellNumber % Looking at every 10th cell
    hold on; 
    plot(cells.rawF(:,i),'k');
    title(['ROI number ', num2str(i)]);
    ylabel('Raw fluorescence'); 
    xlabel('Time (in frames)');
    goodplot;
    pause;
    clf;
end;

%% Set up saving of 'cells' content

%Currently it will simply save the 'cells' variable and its content to a matlab file with the name of the tiff stack in the specified folder as the registered stack was located in
%It might be beneficial to save 'cells' data somewhere else to easier work on personal computer
%Note that the save() function does not work when MIJ has been launched in the same session for some reason, so I use this bc it works. 

uisave('cells',filename); %Will prompt you to confirm where to save the cell variable

clc; clear all; close all;