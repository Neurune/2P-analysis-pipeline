function [imgStack, header] = readTifsRune(dirname)

global HowManyStacks fileNameToGrab; 

imagefiles = dir([dirname '/*.tif']);     

%nfiles = length(imagefiles);

%get image height and width from first file
firstImgPath = [dirname '/' imagefiles(1).name];
firstImg = imread(firstImgPath);
[h w] = size(firstImg);

%determine whether tif is multipage; use scanimage's reader if so.
s = imfinfo(firstImgPath);

imgStack = [];
header = [];

if size(s,1) == 1
    % Assume all files are 1-page tifs
    imgStack = zeros(h,w,nfiles,'uint16');
    for ii=HowManyStacks;
       currentfilename = [dirname '/' imagefiles(ii).name];
       image = imread(currentfilename);
       imgStack(:,:,ii) = image;
    end
else
    
disp('Reading multipage tifs');
    
    % determine total number of images
    nImages = 0;
    for i = HowManyStacks;
        currentfilename = [dirname '/' imagefiles(i).name];
        nPages = size(imfinfo(currentfilename),1);
        nImages = nImages + nPages;
    end
    disp(['Allocating memory for ' num2str(nImages) ' images.']);
    
    % allocate all the memory we need
    imgStack=zeros(h,w,nImages,'uint16');
    
    % read files
    disp(['Reading image data.']);
    ii = 1;
    for i = HowManyStacks
        currentfilename = [dirname '/' imagefiles(i).name];
        nPages = size(imfinfo(currentfilename),1);
        for p = 1:nPages
            imgStack(:,:,ii)=imread(currentfilename,'Index',p);
            ii = ii + 1;
        end
        fileNameToGrab = imagefiles(i).name;
        disp(['Loaded ' imagefiles(i).name]);
    end

end
disp(['Finished loading ' num2str(size(imgStack,3)) ' images.']);


