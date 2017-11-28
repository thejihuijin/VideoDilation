%% 
clear all; close all; clc;
   

%% Load video
% dim_ds = 4; % Downsizing of height & width dimensions
% [rgbVid, fr] = sliceVid('../data/ball_slowmo_crop.mp4', 0, 15, dim_ds);
% 
% slow_speed = 1/8; % If video is in slow-mo change this var
% fr = fr/slow_speed;
% 
% frame_ds = 2; % Downsample the number of frames
% rgbVid = rgbVid(:,:,:,1:2:end);
% fr = fr/frame_ds;
% 
% % Convert to gray & record final dimensions
% [rows, cols, ~, n_frames] = size(rgbVid);
% vid = rgbToGrayVid(rgbVid);
%% Load Video
basefile = 'rs_cat_wall_climb';
filename = strcat(basefile,'.mp4');
dim_ds = 2;

[rgbvid, fr] = sliceVid(filename,0,20,dim_ds);
vid = rgbToGrayVid(rgbvid);
[rows, cols, n_frames] = size(vid);
clear('rgbvid');

%%
resize_vid(filename,'rs_cat_wall_climb.mp4',960,720);
%% Compute Saliency
clc;
iFName = filename;
oFName = strcat(basefile,'_saliency.avi');
oMFName = strcat(basefile,'_saliency.mat');
wsize_s=3; wsize_t=3; wsize_fs=5; wsize_ft=5; 
scaleFactor=1/8; segLength=min(n_frames,150);
detSal(iFName,oFName,oMFName,wsize_s,wsize_t,...
        wsize_fs,wsize_ft,scaleFactor,segLength);
fprintf('Done\n');
clear('regex','wsize_*','segLength','scaleFactor','*FName');
%% Read in Saliency Frames

%saliencyMapHolder = zeros([18,18,448]);
% %saliencyMapTYHolder = zeros([18,18,n_frames]);
% for ii=1:1
%     load(strcat('ws_crop_saliency',num2str(ii),'.mat'));
%     saliencyMapHolder(:,:,(ii-1)*100 + 1:min(ii*100,448)) = salMap(:,:,:);
%     %saliencyMapTYHolder(:,:,(ii-1)*50 + 1:min(ii*50,n_frames)) = salMapTY(:,:,1:2:end);
%     clear('regex','salMap*');
% end
% saliencyMapHolder = saliencyMapHolder(:,:,1:2:end-1);

load(strcat(basefile,'_saliency_1.mat'));
saliencyMapHolder = salMap(:,:,:);

saliencyMapTime = 1/12*(salMapSI + salMapSQ + salMapSY) ...
                + 3/12*(salMapTI + salMapTQ + salMapTY);
%saliencyMapTime = 1/3*(salMapTI + salMapTQ + salMapTY);
saliencyMapS = 1/3*(salMapSI + salMapSQ + salMapSY);
clear('regex','salMap*');
size(saliencyMapHolder)
%% Determine "ground truth"
f=figure; colormap('gray');
for i=1:n_frames
    imagesc(vid(:,:,i))
    %title('Constant framerate')
    title(['Frame ',num2str(i)])
    keydown = waitforbuttonpress;
end
close(f);
%%
gt = ones(1,n_frames)*.5;
gt(80:end) = 1;
save(strcat(basefile,'_gt.mat'),'gt');
% load('dog_and_stuffedDog_gt.mat');
%% Compute optical flow frames
% Horn-Schunck Method
OF = opticalFlowHS();
flow = estimateFlow(OF, vid(:,:,1)); % Use 1st frame as reference
flow_mags = zeros(rows,cols,n_frames);
for i = 1:n_frames
    % I don't know how to initialize a vector of these objects
    flow(i) = estimateFlow(OF, vid(:,:,i));
    
    % Store magnitude frames
    flow_mags(:,:,i) = flow(i).Magnitude;
end


%% Compute dynamic ranges
salminmax = [min(saliencyMapHolder(:)), max(saliencyMapHolder(:))];
ofminmax = [inf,-inf];
ofbounds = zeros(2,n_frames);
for i = 1:n_frames
    temp = min(min(flow(i).Magnitude));
    if temp < ofminmax(1)
        ofminmax(1) = temp;
    end
    ofbounds(1,i) = temp;
    temp = max(max(flow(i).Magnitude));
    if temp > ofminmax(2)
        ofminmax(2) = temp;
    end
    ofbounds(2,i) = temp;
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
    imagesc(tmpSal);%,salminmax);
    title('Saliency Map');
    
    subplot(133);
    imagesc(flow(i).Magnitude,ofminmax);
    title('Optical Flow');
    dur_calc = toc(st);
    pause(1/15 - dur_calc);
end
close(fig);

%% Compute OF Energy
of_total_energy = total_Energy(flow_mags);
of_minkowski = minkowski(flow_mags);
of_min_half = minkowski(flow_mags,.5);
of_five_num = five_num_sum(flow_mags);
of_weight_pool = weight_pool(flow_mags,2);
of_weight_pool_half = weight_pool(flow_mags,.5);

%% Compute Saliency Energy
sal_energy = saliency_EMD_Avg_Energy(saliencyMapHolder);
sal_total_energy = total_Energy(saliencyMapHolder);
sal_minkowski = minkowski(saliencyMapHolder);
sal_min_half = minkowski(saliencyMapHolder,.5);
sal_five_num = five_num_sum(saliencyMapHolder);
sal_weight_pool = weight_pool(saliencyMapHolder,2);
sal_weight_pool_half = weight_pool(saliencyMapHolder,.5);

%% Compute Time-Weighted Saliency Energy
tsal_energy = saliency_EMD_Avg_Energy(saliencyMapTime);
tsal_total_energy = total_Energy(saliencyMapTime);
tsal_minkowski = minkowski(saliencyMapTime);
tsal_min_half = minkowski(saliencyMapTime,.5);
tsal_five_num = five_num_sum(saliencyMapTime);
tsal_weight_pool = weight_pool(saliencyMapTime,2);
tsal_weight_pool_half = weight_pool(saliencyMapTime,.5);




%% Compare
figure; hold on
if exist('gt','var')
    plot(gt,'DisplayName','Ground Truth');
end
plot(of_total_energy,'DisplayName','OF Total');
plot(of_minkowski,'DisplayName','OF Minkowski');
plot(of_min_half,'DisplayName','OF Minkowski 1/2');
plot(of_five_num,'DisplayName','OF 5 Num Sum');
plot(of_weight_pool,'DisplayName','OF Weight Pool 2');
plot(of_weight_pool_half,'DisplayName','OF Weight Pool 1/2');
legend();

title('OF Energy Functions');
figure; hold on
if exist('gt','var')
    plot(gt,'DisplayName','Ground Truth');
end
plot(sal_energy,'DisplayName','sal EMD');
plot(sal_total_energy,'DisplayName','sal Total');
plot(sal_minkowski,'DisplayName','sal Minkowski');
plot(sal_min_half,'DisplayName','sal Minkowski 1/2');
plot(sal_five_num,'DisplayName','sal 5 Num Sum');
plot(sal_weight_pool,'DisplayName','sal Weight Pool 2');
plot(sal_weight_pool_half,'DisplayName','sal Weight Pool 1/2');
legend();
title('Sal Energy Functions');

figure; hold on
if exist('gt','var')
    plot(gt,'DisplayName','Ground Truth');
end
plot(tsal_energy,'DisplayName','tsal EMD');
plot(tsal_total_energy,'DisplayName','tsal Total');
plot(tsal_minkowski,'DisplayName','tsal Minkowski');
plot(tsal_min_half,'DisplayName','tsal Minkowski 1/2');
plot(tsal_five_num,'DisplayName','tsal 5 Num Sum');
plot(tsal_weight_pool,'DisplayName','tsal Weight Pool 2');
plot(tsal_weight_pool_half,'DisplayName','tsal Weight Pool 1/2');
legend();
title('Time Saliency Energy Functions');

%% Compute Saliency Masking
masked_flow = zeros(size(flow_mags));
masked_tsal = zeros(size(saliencyMapTime));
for i=1:n_frames
    sal_mask = imresize(saliencyMapHolder(:,:,i),[rows, cols],'bilinear');
    sal_mask = imbinarize(sal_mask);
    
    masked_flow(:,:,i) = flow_mags(:,:,i).*sal_mask;
    
    sal_mask = imbinarize(saliencyMapHolder(:,:,i));
    masked_tsal(:,:,i) = saliencyMapTime(:,:,i).*sal_mask;
end
%% Compute Masked Time-Weighted Saliency Energy
mtsal_energy = saliency_EMD_Avg_Energy(masked_tsal);
mtsal_total_energy = total_Energy(masked_tsal);
mtsal_minkowski = minkowski(masked_tsal);
mtsal_min_half = minkowski(masked_tsal,.5);
mtsal_five_num = five_num_sum(masked_tsal);
mtsal_weight_pool = weight_pool(masked_tsal,2);
mtsal_weight_pool_half = weight_pool(masked_tsal,.5);
%% Compute Masked OF Energy
mof_energy = saliency_EMD_Avg_Energy(masked_flow);
mof_total_energy = total_Energy(masked_flow);
mof_minkowski = minkowski(masked_flow);
mof_min_half = minkowski(masked_flow,.5);
mof_five_num = five_num_sum(masked_flow);
mof_weight_pool = weight_pool(masked_flow,2);
mof_weight_pool_half = weight_pool(masked_flow,.5);
%%
sal_masks = zeros(size(vid));
sal_masksS = zeros(size(vid));
sal_masksT = zeros(size(vid));
for i=1:n_frames
    sal_masks(:,:,i) = imbinarize(imresize(saliencyMapHolder(:,:,i),[rows,cols],'bilinear')).*vid(:,:,i);
    sal_masksS(:,:,i) = imbinarize(imresize(saliencyMapS(:,:,i),[rows,cols],'bilinear')).*vid(:,:,i);
    sal_masksT(:,:,i) = imbinarize(imresize(saliencyMapTime(:,:,i),[rows,cols],'bilinear')).*vid(:,:,i);

end
%%
figure; colormap('gray');
for i=1:n_frames
    subplot(121)
    imagesc(sal_masksT(:,:,i));
    title('Overall Saliency');
    
    subplot(122);
    imagesc(sal_masksS(:,:,i));
    title('Spatial Saliency');
    pause(1/10);
end
%%
figure; hold on
if exist('gt','var')
    plot(gt,'DisplayName','Ground Truth');
end
plot(mtsal_total_energy,'DisplayName','mtsal Total');
plot(mtsal_minkowski,'DisplayName','mtsal Minkowski');
plot(mtsal_min_half,'DisplayName','mtsal Minkowski 1/2');
plot(mtsal_five_num,'DisplayName','mtsal 5 Num Sum');
plot(mtsal_weight_pool,'DisplayName','mtsal Weight Pool 2');
plot(mtsal_weight_pool_half,'DisplayName','mtsal Weight Pool 1/2');
legend();
title('Masked Time Saliency Energy Functions');

figure; hold on
if exist('gt','var')
    plot(gt,'DisplayName','Ground Truth');
end
plot(mof_energy,'DisplayName','mof EMD');
plot(mof_total_energy,'DisplayName','mof Total');
plot(mof_minkowski,'DisplayName','mof Minkowski');
plot(mof_min_half,'DisplayName','mof Minkowski 1/2');
plot(mof_five_num,'DisplayName','mof 5 Num Sum');
plot(mof_weight_pool,'DisplayName','mof Weight Pool 2');
plot(mof_weight_pool_half,'DisplayName','mof Weight Pool 1/2');
legend();
title('Masked Optical Flow Energy Functions');
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
%%
energy_normal = 1-weight_pool(masked_flow,.5);

%% Scale framerate to desired range
% Framerate normalized to range between min_fr and min_fr + range_fr
min_fr = 15;
range_fr = 60;
fr_scaled = energy_normal*range_fr + min_fr;
fr_scaled2 = 30.*2.^((-1+2*energy_normal));

figure; hold on
plot(fr_scaled,'DisplayName','Linear');
plot(fr_scaled2,'DisplayName','Exp');
legend()

% Alternatively, scale so that the mean framerate = original framerate
% energy_normal2 = energy_normal / mean(energy_normal);
% fr_scaled2 = min_fr + (fr - min_fr)*energy_normal2;

%% Smooth, renormalize, and scale the energy-normalized framerate
mov_avg_window = 15;

%energy_smooth = movmean(energy_normal, mov_avg_window);
%energy_smooth = (energy_smooth - min(energy_smooth))/range(energy_smooth);
%fr_smoothed = energy_smooth*range_fr + min_fr;
fr_smoothed = movmean(fr_scaled, mov_avg_window);
fr_smoothed2 = movmean(fr_scaled2, mov_avg_window);
% Alternatively, scale so that the mean framerate = original framerate
% energy_smoothed2 = movmean( energy_normal, mov_avg_window);
% energy_smoothed2 = energy_smoothed2 / mean(energy_smoothed2);
% fr_smoothed2 = min_fr + (fr - min_fr)*energy_smoothed2;

%% Quantize framerate, adjust framerate changes, re-smooth framerate
q_levels = 5;
fr_quant = quantFR( fr_smoothed, q_levels );
fr_q_adj = adjustFR( fr_quant, 0.1, fr );

% Smooth the adjusted framerate
mov_avg_window = 5;
fr_q_adj_smooth = movmean( fr_q_adj, mov_avg_window );
fr_q_adj_smooth = movmean( fr_q_adj_smooth, mov_avg_window );

figure
plot(1:n_frames,fr_smoothed, 1:n_frames,fr_q_adj, 1:n_frames, fr_q_adj_smooth)
legend('Smoothed', 'Adjusted, Quantized', 'Adjusted, Quantized, Smoothed')

%% Compute playback vector from framerate vector
scaled_playback = fr2playback(fr_scaled, fr);
smoothed_playback = fr2playback(fr_smoothed, fr);
adjusted_smooth_playback = fr2playback(fr_q_adj_smooth, fr);
%% Play Video
figure
playDilatedFrames( vid, adjusted_smooth_playback, fr, fr_q_adj_smooth )
%% Play Videos
figure, colormap('Gray');
currTime = 0;
for i = 1:n_frames
    imagesc(vid(:,:,i))
    title(['Original - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/30 );
    currTime = currTime + 1/30;
end

%%
figure, colormap('Gray')
currTime = 0;
for i = 1:n_frames
    imagesc(vid(:,:,i))
    title(['Optical Flow @ ',sprintf('%i',round(of_fr(i))),' fps - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/of_fr(i) );
    currTime = currTime + 1/of_fr(i);
end
%%
figure, colormap('Gray')
currTime = 0;
for i = 1:n_frames
    imagesc(vid(:,:,i))
    title(['Saliency @ ',sprintf('%i',round(sal_fr(i))),' fps - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/sal_fr(i) );
    currTime = currTime + 1/sal_fr(i);
end





























