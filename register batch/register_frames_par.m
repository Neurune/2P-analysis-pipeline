function [shift_x,shift_y]=register_frames_par(data, template, method, usfac)
% aligns images in stack to a template using xcorr 
% input :
%          data : the stack of images ; a matrix of the form w x h x n,
%                  where w: width, h: height and n, nb of images in stack.
%
%          template : the image used as reference for the registration
%                     if no template is provided, align to the mean
%          
%          method : a string that can be 'fft'(original method using xcorr)
%                   or 'dft' (faster, but can also go to subpixel accuracy)
%
%          usfac : upsampling factor: determines the shift to the precision
%                  ( pixelsize/usfac ), per default usfac = 1 
%                  'fft' method works implicitly with usfac = 1 
%
% output :
%         shift = [shift_x,shift_y] vector of size n x 2 of 2D-shifts 
%                 expressed in pixel units
%
% -----------------------------------12/04/2012, Raphael Thierry, FMI/FAIM

% default input values :
if nargin<2, template=mean(data,3); end
if nargin<3, method='fft';          end
if nargin<4, usfac = 1;             end

size_data = size(data);

% user-defined values (hard-coded)
%------------------------------------
% chunks for parallel processing
n_chunks = 1; % data is separated in n_chunks ( original value = 20;)
              % place here the number of core of the machine (4, 8,...)
              
              %RR: Can this work with > 1 n_chunks?
              
chunk_size = round(size_data(3)/n_chunks);

% refreshs the template with mean(nbim_temp images), every nbim_refresh im
nbim_temp = 5; % template = mean over nbim_temp number of images 
nbim_refresh = round(chunk_size/10); % pace of template refreshment

perc_boundary = 0.05; % percentage of raws/columns cropped on im boundaries; play around with this value!!
template=template-mean2(template);
% determine how much of the images to use for alignment, the larger the
% boundary the less pixels are used for alignment and the faster the
% algorithm runs.
boundary=round(perc_boundary*max(size(data(:,:,1))));
template=template(boundary+1:end-boundary,boundary+1:end-boundary);

chunk_inds=[1:chunk_size:(n_chunks-1)*chunk_size+1 size(data,3)+1];


% enabling parallel processing
%------------------------------
% if  matlabpool('size')==0
%     try
%         matlabpool open
%     end
% end
%------------------------------

if strcmp(method,'fft')

    % output variable allocation
    shift_x=zeros(size(data,3),1);
    shift_y=zeros(size(data,3),1);
    %c=fft_xcorr(template,data(boundary+1:end-boundary,boundary+1:end-boundary,ind));

    F_template_size = [size(template,1) size(template,2)];
    F_in_size = [size(data,1)-2*boundary size(data,2)-2*boundary];
    outsize = F_template_size + F_in_size - 1;

    low_pass_thresh=round(min(F_template_size)/4);
    high_pass_thresh=round(min(F_template_size)/40);

    template = im2double(template);

    % fourier transform and band pass filter the template
    F_template = fft2(rot90(template,2),outsize(1),outsize(2));
    F_template([1:high_pass_thresh size(F_template,1)-high_pass_thresh+2:size(F_template,1)],:)=0;
    F_template(:,[1:high_pass_thresh size(F_template,2)-high_pass_thresh+2:size(F_template,2)])=0;
    F_template(low_pass_thresh+2:end-low_pass_thresh,:)=0;
    F_template(:,low_pass_thresh+2:end-low_pass_thresh)=0;

    for knd = 1:n_chunks
        
        disp(['Starting Chunk ' num2str(knd) ' of ' num2str(n_chunks) ])
        data_chunk=data(:,:,chunk_inds(knd):chunk_inds(knd+1)-1);
        shift_x_chunk=zeros(n_chunks,size(data_chunk,3),1);
        shift_y_chunk=zeros(size(data_chunk,3),1);
        
        tic;
        parfor ind=1:size(data_chunk,3) 
            
            sc_in=double(data_chunk(:,:,ind))-mean(mean(data_chunk(:,:,ind)));
            % fourier transform and band pass filter current frame
            F_in = fft2(sc_in(boundary+1:end-boundary,boundary+1:end-boundary),outsize(1),outsize(2));
            F_in(low_pass_thresh+2:end-low_pass_thresh,:)=0;
            F_in(:,low_pass_thresh+2:end-low_pass_thresh)=0;
            F_in([1:high_pass_thresh size(F_in,1)-high_pass_thresh+2:size(F_in,1)],:)=0;
            F_in(:,[1:high_pass_thresh size(F_in,2)-high_pass_thresh+2:size(F_in,2)])=0;

            % calculate the cross correlation 
            c = real(ifft2(F_template .* F_in));
            c=c(boundary+1:end-boundary,boundary+1:end-boundary);

            % determine the indices of the maximum
            [not_used,max_ind]=max(c(:));
            [max_x,max_y] = ind2sub(size(c),max_ind);

            % calculate the x and y shifts.
            shift_x_chunk(ind)=max_x-ceil(size(c,1)/2);
            shift_y_chunk(ind)=max_y-ceil(size(c,2)/2);
        end
        toc %disp(['Chunk ' num2str(knd) ' of ' num2str(n_chunks) ' - Elapsed time: ' num2str(round(toc/60)) ' min'])
        shift_x(chunk_inds(knd):chunk_inds(knd+1)-1)=shift_x_chunk;
        shift_y(chunk_inds(knd):chunk_inds(knd+1)-1)=shift_y_chunk;

    end

elseif strcmp(method,'dft')
    
     %fft2(mean(im2double(data_chunk(:,:,1:nbim_temp)),3));
    % output variable allocation
    shift_x=zeros(n_chunks,chunk_size);
    shift_y=zeros(n_chunks,chunk_size); 
    data = reshape(data(boundary+1:end-boundary, boundary+1:end-boundary, :),[size_data(1)-2*boundary, size_data(2)-2*boundary, chunk_size, n_chunks]);
    
    tic
    
    parfor knd = 1:n_chunks
        
        fft2_template = fft2(template);
        %disp(['Starting Chunk ' num2str(knd) ' of ' num2str(n_chunks) ])
        %data_chunk = data(boundary+1:end-boundary, boundary+1:end-boundary, chunk_inds(knd):chunk_inds(knd)+(chunk_size-1));
        data_chunk = squeeze(data(:,:,:,knd));
        shifts=zeros(chunk_size,4);
        
        for ind = 1:chunk_size 
                      
            if ~mod(ind,nbim_refresh) 
                % refresh the template by aligning the last nbim_temp imag
                temp_ali = shift_data(double(data_chunk(:,:,ind-nbim_temp:ind-1)),shifts(ind-nbim_temp:ind-1,3),shifts(ind-nbim_temp:ind-1,4));
                % ...and take the fft from the mean
                fft2_template = fft2(im2double(mean(temp_ali,3)));     
            end
            
            % calculates registration via DFT :
             shifts(ind,:)= dftregistration(fft2(double(data_chunk(:,:,ind))), fft2_template, usfac);
                   
        end
        tmp = shifts(:,3);
        shift_x(knd,:)= tmp'; 
        tmp = shifts(:,4);
        shift_y(knd,:)= tmp';
            
    end
    
    toc
    shift_x = shift_x';
    shift_x = shift_x(:);
    shift_y = shift_y';
    shift_y = shift_y(:);
    
end
