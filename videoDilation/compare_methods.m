%% 
clear all; close all; clc;

% vid = rgbToGrayVid(rgbVid);
%% Load Video
basefile = 'rs_cat_dunk';
filename = strcat(basefile,'.mp4');
dim_ds = 4;

[rgbvid, fr] = sliceVid(filename,0,20,dim_ds);
vid = rgbToGrayVid(rgbvid);
[rows, cols, n_frames] = size(vid);
clear('rgbvid');

%%
resize_vid('cat_dunk.mp4','rs_cat_dunk.mp4',720,1200);
%% Compute Saliency
clc;
iFName = filename;
oFName = strcat(basefile,'_saliency.avi');
oMFName = strcat(basefile,'_saliency.mat');
wsize_s=3; wsize_t=3; wsize_fs=5; wsize_ft=5; 
scaleFactor=1/8; segLength=min(n_frames,50);
%detSal(iFName,oFName,oMFName,wsize_s,wsize_t,...
%        wsize_fs,wsize_ft,scaleFactor,segLength);
[smap, smapt] = compute_saliency(iFName,wsize_s,wsize_t,wsize_fs,wsize_ft,scaleFactor,segLength);
fprintf('Done\n');
clear('regex','wsize_*','segLength','scaleFactor','*FName');

size(smap)
size(smapt)
%%
fig = figure;
for i = 1:n_frames
    subplot(121);
    imagesc(smap(:,:,i))
    %title('Constant framerate')
    title(['Original @ ' num2str(15) ' fps - ' ...
        sprintf('%.2f',i/15) ' seconds elapsed'])
    
    subplot(122);
    %tmpSal=imresize(saliencyMapTime(:,:,i),[rows, cols],'bilinear');
    imagesc(smapt(:,:,i));%,salminmax);
    title('Saliency Map');
    pause(1/15);
end
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
% f=figure; colormap('gray');
% for i=1:n_frames
%     imagesc(vid(:,:,i))
%     %title('Constant framerate')
%     title(['Frame ',num2str(i)])
%     keydown = waitforbuttonpress;
% end
% close(f);
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

%% Compare
% figure; hold on
% if exist('gt','var')
%     plot(gt,'DisplayName','Ground Truth');
% end
% plot(of_total_energy,'DisplayName','OF Total');
% plot(of_minkowski,'DisplayName','OF Minkowski');
% plot(of_min_half,'DisplayName','OF Minkowski 1/2');
% plot(of_five_num,'DisplayName','OF 5 Num Sum');
% plot(of_weight_pool,'DisplayName','OF Weight Pool 2');
% plot(of_weight_pool_half,'DisplayName','OF Weight Pool 1/2');
% legend();


%% Compute Saliency Masking
masked_flow = zeros(size(flow_mags));
for i=1:n_frames
    sal_mask = imresize(saliencyMapHolder(:,:,i),[rows, cols],'bilinear');
    sal_mask = imbinarize(sal_mask);
    
    masked_flow(:,:,i) = flow_mags(:,:,i).*sal_mask;
end

clear('sal_mask');


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
%%
energy_normal = 1-weight_pool(masked_flow,.5);

%% Scale framerate to speed up
% Adjust framerate to exponential
scale = 1;
fr_scaled = fr.*2.^(scale*(-1+2*energy_normal));

% Smooth Adjusted Frame Rate
mov_avg_window = 15;
fr_smoothed = movmean(fr_scaled, mov_avg_window);


%% Quantize framerate, adjust framerate changes, re-smooth framerate
q_levels = 5;
%fr_quant = quantFR( fr_smoothed, q_levels );
fr_q_adj = adjustFR( fr_smoothed, 0.2, fr );

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


























