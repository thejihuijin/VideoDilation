%% 
clear all; close all; clc;

%% Resize video for saliency
% dimensions need to be divisible by 1/scaleFactor*wsize_s*wsize_fs
% for current parameters, this value is 120
% resize_vid('cat_dunk.mp4','rs_cat_dunk.mp4',720,1200);
%% Load Video
basefile = 'rs_cat_dunk';
filename = strcat(basefile,'.mp4');
dim_ds = 4;

[rgbvid, fr] = sliceVid(filename,0,20,dim_ds);
vid = rgbToGrayVid(rgbvid);
[rows, cols, n_frames] = size(vid);
clear('rgbvid','dim_ds');


%% Compute Saliency
clc;
wsize_s=3; wsize_t=3; wsize_fs=5; wsize_ft=5; 
scaleFactor=1/8; segLength=100;

[saliencyMapHolder, saliencyMapTime] = compute_saliency(filename,wsize_s,wsize_t,wsize_fs,wsize_ft,scaleFactor,segLength);
fprintf('Done\n');

clear('regex','wsize_*','segLength','scaleFactor','*FName');

%% Compute optical flow frames
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
fig = figure;
for i = 1:n_frames
    st = tic;
    
    subplot(131);
    imagesc(vid(:,:,i))
    %title('Constant framerate')
    title(['Original @ ' num2str(15) ' fps - ' ...
        sprintf('%.2f',i/15) ' seconds elapsed'])
    
    subplot(132);
    tmpSal=imresize(saliencyMapTime(:,:,i),[rows, cols],'bilinear');
    imagesc(tmpSal, salminmax);
    title('Saliency Map');
    
    subplot(133);
    imagesc(flow_mags(:,:,i),ofminmax);
    title('Optical Flow');
    dur_calc = toc(st);
    pause(1/15 - dur_calc);
end
close(fig);

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
close all;
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
%%
clc; close all;
%% Compute Frame Rates from energy
% select energy function
energy = compute_energy(flow_mags,saliencyMapHolder,saliencyMapTime,'MOF','WP');

% convert to fr
time_padded_fr=energy2fr(energy,fr,.2);

%% Compute playback vector from framerate vector
% resample frames based on new frame rate
adjusted_smooth_playback = fr2playback(time_padded_fr, fr);
%% Play Video
figure
playDilatedFrames( vid, adjusted_smooth_playback, fr, time_padded_fr )


























