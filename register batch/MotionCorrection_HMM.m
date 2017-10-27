function [dl,dk,max_num_of_displacements]=MotionCorrection_HMM(I,R)
% this function is an implementation of the dombeck&tank HMM alignment
% algorithm
% I: raw image (to be aligned)
% R: template

I=double(I)';
R=R';

%filter 
% mI=mean(I(:));
% mR=mean(R(:));

% high_pass_thresh=round(min(size(R))/80);
% 
% I = fft2(I);
% I([2:high_pass_thresh size(I,1)-high_pass_thresh+2:size(I,1)],:)=0;
% I(:,[2:high_pass_thresh size(I,2)-high_pass_thresh+2:size(I,2)])=0;
% I = abs(ifft2(I));
% 
% R = fft2(R);
% R([2:high_pass_thresh size(R,1)-high_pass_thresh+2:size(R,1)],:)=0;
% R(:,[2:high_pass_thresh size(R,2)-high_pass_thresh+2:size(R,2)])=0;
% R = abs(ifft2(R));


% parameters
max_num_of_displacements=5;
displacement_step=1; % in pixels

% increase this value to increase the likelihood of larger jumps between
% successive scan lines
trans_prob_space_const=2*1e5;

% do not change this value
gamma=6.12e-006*trans_prob_space_const^2;
frame_ind=1;
number_of_frames=1;
color_scale_boundary=1;

% initialize variables
transition_probability=[];
image_width=size(I,1);
image_height=size(I,2);

% define the transition probabilities
for ind=1:4*max_num_of_displacements+1
    for knd=1:4*max_num_of_displacements+1
        a(ind,knd)=1/(2*pi*trans_prob_space_const^2)*exp(-sqrt((ind-1-2*max_num_of_displacements)^2+(knd-1-2*max_num_of_displacements)^2)/trans_prob_space_const*displacement_step);
    end
end




% initialize the best path probability
j=[max_num_of_displacements*displacement_step+1 : image_width-max_num_of_displacements*displacement_step];
k=max_num_of_displacements*displacement_step+1;

for ind=1:2*max_num_of_displacements+1
    for knd=1:2*max_num_of_displacements+1
        dx=displacement_step*(ind-1-max_num_of_displacements);
        dy=displacement_step*(knd-1-max_num_of_displacements);
        V(ind,knd,1)=a(ind+max_num_of_displacements,knd+max_num_of_displacements)*gamma*sum(I(j,k,frame_ind).*log(R(j+dx,k+dy))-R(j+dx,k+dy));
        %V(ind,knd,1)=gamma*sum(I(j,k,frame_ind).*log(R(j+dx,k+dy))-R(j+dx,k+dy));
        B(ind,knd,1,[1:2])=0;
    end
end

a(:,[1:2*max_num_of_displacements-1])=0;
a(:,[2*max_num_of_displacements+3:end])=0;
a([1:2*max_num_of_displacements-1],:)=0;
a([2*max_num_of_displacements+3:end],:)=0;

% hf=figure;
% set(hf,'menubar','none','position',[600 600 300 50])
% ha=area([0 0],[1 1],'facecolor','r','edgecolor','none');
% xlim([0 1]);
% axis off
% title('Calculating...')
% drawnow;

% forward algorithm
for line_ind=2:image_height-2*max_num_of_displacements*displacement_step
    k=max_num_of_displacements*displacement_step+line_ind;
    for ind=1:1:2*max_num_of_displacements+1
        for knd=1:2*max_num_of_displacements+1
            dx=displacement_step*(ind-1-max_num_of_displacements);
            dy=displacement_step*(knd-1-max_num_of_displacements);
            trans_prob_tmp=circshift(a,[ind, knd]-max_num_of_displacements-1);
            %keyboard
            trans_prob_tmp=trans_prob_tmp(max_num_of_displacements+1:end-max_num_of_displacements,max_num_of_displacements+1:end-max_num_of_displacements);
            [tmp_V, tmp_B] = max( V(:,:,line_ind-1).*trans_prob_tmp*gamma*sum(I(j,k,frame_ind).*log(R(j+dx,k+dy))-R(j+dx,k+dy)) );
            [V(ind,knd,line_ind) B(ind,knd,line_ind,2)] = max(tmp_V);
            B(ind,knd,line_ind,1)=tmp_B(B(ind,knd,line_ind,2));
        end
    end
    %     figure(1);imagesc(V(:,:,line_ind))
    %     keyboard
    
%     set(ha,'XData',[0 line_ind/(image_height-2*max_num_of_displacements*displacement_step)])
%     xlim([0 1]);
%     drawnow
end
%close(hf)

% backtracking the hidden state sequence

[tmp_V, tmp_delta]=max(V(:,:,line_ind));
[max_V, dl(line_ind)]=max(tmp_V);
dk(line_ind)=tmp_delta(dl(line_ind));

for line_ind=image_height-2*max_num_of_displacements*displacement_step-1:-1:1
    dk(line_ind)=B(dk(line_ind+1),dl(line_ind+1),line_ind+1,1);
    dl(line_ind)=B(dk(line_ind+1),dl(line_ind+1),line_ind+1,2);
end

dk=dk-max_num_of_displacements-1;
dl=dl-max_num_of_displacements-1;
















