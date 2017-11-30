%% 
clear all; close all; clc;
%% Define parameters
wsize_s=3; wsize_t=3; wsize_fs=5; wsize_ft=5; 
scaleFactor=1/8; segLength=100;

dim_ds = 2;
%% Check Video size for saliency algorithm
% NOTE: Path must be absolute or correct relative path to video (i.e. you
% must be in the correct directory for relative path to work).
% If the video dimension does not meet the criteria of the saliency 
% algorithm, a new video will be generated
filename = 'data/Reddit_Videos/zaboomafoo.mp4';
filename = check_video(filename, wsize_s*wsize_fs/scaleFactor);

%% Load Video
[rgbvid, fr] = sliceVid(filename,0,20,dim_ds);
vid = rgbToGrayVid(rgbvid);
[rows, cols, n_frames] = size(vid);
clear('rgbvid','dim_ds');


%% Compute Saliency
fprintf('Computing Saliency\n');
[saliencyMapHolder, saliencyMapTime] = compute_saliency(filename,wsize_s,wsize_t,wsize_fs,wsize_ft,scaleFactor,segLength);
fprintf('Done\n');

saliencyMapHolder = saliencyMapHolder(:,:,1:n_frames);
saliencyMapTime = saliencyMapTime(:,:,1:n_frames);
%% Compute optical flow frames
fprintf('Computing Optical Flow\n');
flow_mags = compute_OF(vid);


%% Compute dynamic ranges
% For display purposes only
salminmax = [min(saliencyMapHolder(:)), max(saliencyMapHolder(:))];
ofminmax = [inf,-inf];
for i = 1:n_frames
    temp = min(min(flow_mags(:,:,i)));
    if temp < ofminmax(1)
        ofminmax(1) = temp;
    end
    temp = max(max(flow_mags(:,:,i)));
    if temp > ofminmax(2)
        ofminmax(2) = temp;
    end
end
clear('temp');
%% Display Saliency and Optical Flow
fprintf('Displaying Original Video\n');
figure; colormap gray;
for i = 1:n_frames
    st = tic;
    
    %subplot(131);
    imagesc(vid(:,:,i))
    %title('Constant framerate')
    title(['Original @ ' num2str(fr) ' fps - ' ...
        sprintf('%.2f',i/fr) ' seconds elapsed'])
    
%     subplot(132);
%     tmpSal=imresize(saliencyMapTime(:,:,i),[rows, cols],'bilinear');
%     imagesc(tmpSal, salminmax);
%     title('Saliency Map');
%     
%     subplot(133);
%     imagesc(flow_mags(:,:,i),ofminmax);
%     title('Optical Flow');
    dur_calc = toc(st);
    pause(1/fr - dur_calc);
end


%%
%
%
%                 Compare Energy Functions
%
%
%% Compute Energy Functions
% Optical Flow
of_minkowski = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'OF','MINK');
of_five_num = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'OF','FNS');
of_weight_pool_half = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'OF','WP');

% Time Weighted Saliency
tsal_minkowski = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'TSAL','MINK');
tsal_five_num = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'TSAL','FNS');
tsal_weight_pool_half = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'TSAL','WP');

% Masked OF
mof_minkowski = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'MOF','MINK');
mof_five_num = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'MOF','FNS');
mof_weight_pool_half = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'MOF','WP');
%% Final Graph

xvals = 1:3:n_frames;
figure; hold on
if exist('gt','var')
    plot(gt,'DisplayName','Ground Truth');
end
% Plot OF
plot(xvals,of_minkowski(1:3:end),'--r*','DisplayName','OF Minkowski');
plot(xvals,of_weight_pool_half(1:3:end),'--ro','DisplayName','OF Weight Pool 1/2');
plot(xvals,of_five_num(1:3:end),'--rx','DisplayName','OF Five Num Sum');

% Plot Time Saliency
plot(xvals,tsal_minkowski(1:3:end),'--g*','DisplayName','TS Minkowski');
plot(xvals,tsal_weight_pool_half(1:3:end),'--go','DisplayName','TS Weight Pool 1/2');
plot(xvals,tsal_five_num(1:3:end),'--gx','DisplayName','TS Five Num Sum');


% Plot Masked OF
plot(xvals,mof_minkowski(1:3:end),'--b*','DisplayName','MOF Minkowski');
plot(xvals,mof_weight_pool_half(1:3:end),'--bo','DisplayName','MOF Weight Pool 1/2');
plot(xvals,mof_five_num(1:3:end),'--bx','DisplayName','MOF Five Num Sum');

legend()
title('Comparison of Various Energy Functions')
xlabel('Frame Number')
%%
clear('regex','of_*','sal_*','tsal_*','mof_*','mtsal_*');

%%
%
%
%                 PLAYBACK VIDEOS
%
%

%% Compute Frame Rates from energy
% select energy function
energy = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'MOF','WP');

% convert to fr
time_padded_fr=energy2fr(energy,fr,.2);

%% Compute playback vector from framerate vector
% resample frames based on new frame rate
adjusted_smooth_playback = fr2playback(time_padded_fr, fr);
%% Display Energy and Framerates
figure;
subplot(121)
plot(1-energy);
xlabel('Frame Number');
ylabel('Inverted Energy');
title('Inverted Energy Graph of Video');

subplot(122);
plot(time_padded_fr);
xlabel('Frame Number');
ylabel('Frame Rate');
title('Time Dilated Frame Rate');
%% Play Video
figure
playDilatedFrames( vid, adjusted_smooth_playback, fr, time_padded_fr )


























