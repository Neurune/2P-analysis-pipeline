%% Script for registering 2P image stack versin 2.0.

%Registration method is based on advice from Dan Wilson (MPFI) and uses the NoRMCorre method:
%Documentation: https://github.com/simonsfoundation/NoRMCorre

%Make sure to use the Bio-format GUI for saving a non-compressed, unsigned
%tif stack of the imaging data after loading in the raw stack, otherwise it will not load into Matlab.

%Comments: 

%Will save the tif stacks with sampling rate of 15 Hz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all; clc;

%% Get the dir for the unregistered 2P image stack and set output dir
inputDir = '/Users/RuneRasmussen/Desktop/inputDir/'; %Input dir
outputDir = '/Users/RuneRasmussen/Desktop/outputDir/'; %Output dir
cd(inputDir); %Jump to this directory

%% How many tif stacks in inputDir and file names
imagefiles = dir([inputDir '/*.tif']);     
nfiles = length(imagefiles);
Method = 1;
downfactor = 2; %i.e. from 30 to 15 Hz - this is what Hofer and Wilson did
time_stated = clock;

global HowManyStacks fileNameToGrab;
for HowManyStacks = 1:length(imagefiles);
    
% Load in the stack(s)
warning off;
disp('reading in data');
data = readTifsRune(inputDir); %Loading stacks into data array, might take long for large stacks, but it will load
filename{HowManyStacks} = fileNameToGrab;
warning on;

% Method = 1; %1 = Keller method, 2 = NoRMCorre method;
 
if Method == 2;

% Perform motion registration
data = single(data); %Convert to single precision 
T = size(data,ndims(data)); %Get dimensions. You might need to change preferences in Matlab if stack is very large

warning off;
% now try non-rigid motion correction (uses parallel pool) - This takes a long tim!!!
% needs to optimize parameters for getting better registration; Ask Dan Wilson for good parameters
options_nonrigid = NoRMCorreSetParms('d1',size(data,1),'d2',size(data,2),'grid_size',[32,32],'mot_uf',4,'bin_width',50,'max_shift',15,'max_dev',3,'us_fac',1,'init_batch',200, 'shifts_method', 'FFT','iter',3);
[M2,shifts2,template2] = normcorre_batch(data,options_nonrigid);
warning on;

%% Compute metrics for how well the registration performed
nnY = quantile(data(:),0.005);
mmY = quantile(data(:),0.995);
T = length(cY);
    
%% plot a movie with the results

%Show downsampled movie - i.e. from 30 Hz to 6 Hz 
% data_ds = downsample_data(data, 'time',downfactor,1);
% M2_ds = downsample_data(M2, 'time',downfactor,1);
% 
% figure;
% for t = 1:1:(T/downfactor);
%     subplot(121);imagesc(data_ds(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     subplot(122);imagesc(M2_ds(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
%     drawnow;
%     pause(0.02);
% end

else
    
%Determine the template to use for determining XY shifts
template = uint16(mean(data(:,:,1:200),3)); %Currently using first 200 frames as template     

%Get the XY shifts
[shift_x,shift_y] = register_frames_par(data, template, 'fft',1); %

%Shift input data using designated value in X and Y dimension 
[data_registered] = shift_data(data,shift_x,shift_y);
data_registered = downsample_data(data_registered, 'time', downfactor,1);

% T = size(data,ndims(data));
% nnY = quantile(data_ds(:),0.005);
% mmY = quantile(data_ds(:),0.995);

% figure;
% for t = 1:1:(T/downfactor);
%     subplot(121);imagesc(data_ds(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,(T/downfactor)),'fontweight','bold','fontsize',14); colormap('bone')
%     subplot(122);imagesc(registered_ds(:,:,t),[nnY,mmY]); xlabel('Keller Method','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,(T/downfactor)),'fontweight','bold','fontsize',14); colormap('bone')
%     drawnow;
%     pause(0.02);
% end
    
end

%%% Saving the registered stack
fixedEnding = '_registered15Hz.tif'; %make the X Hz dynamic based on downsample factor
filenameToSave = filename{HowManyStacks}; filenameToSave = filenameToSave(1:end-4);
outputDirFilename = strcat(outputDir,filenameToSave,fixedEnding);
 
% This can only write up to about 8500 frames. If larger stacks need to be saved we should split it into saved chunks
[x y z] = size(data_registered);
disp(['I am now writing ' filenameToSave 'to ' outputDir]);
for i=1:z
    imwrite(uint16(squeeze(data_registered(:,:,i))),outputDirFilename,'tif','writemode','append');
end
disp(['I am  done writing ' filenameToSave 'to' outputDir]);

%clearvars data data_registered; %To clear up space in Workspace

end
disp(['I am done registering your 2P imaging stacks']);
time_finished = clock;
