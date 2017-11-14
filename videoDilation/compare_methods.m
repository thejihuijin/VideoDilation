%% 
clear all; close all; clc;
   
%% Compute Saliency
% 
% iFName = 'ball_slowmo_crop.mp4';
% oFName = 'ball_slowmo_crop_saliency.mpg';
% oMFName = 'ball_slowmo_crop_saliency.mat';
% wsize_s=3; wsize_t=3; wsize_fs=5; wsize_ft=5; 
% scaleFactor=1/8; segLength=100;
% detSal(iFName,oFName,oMFName,wsize_s,wsize_t,...
%         wsize_fs,wsize_ft,scaleFactor,segLength);
%% Load video
dim_ds = 4; % Downsizing of height & width dimensions
[rgbVid, fr] = sliceVid('../data/ball_slowmo_crop.mp4', 0, 15, dim_ds);

slow_speed = 1/8; % If video is in slow-mo change this var
fr = fr/slow_speed;

frame_ds = 2; % Downsample the number of frames
rgbVid = rgbVid(:,:,:,1:2:end);
fr = fr/frame_ds;

% Convert to gray & record final dimensions
[rows, cols, ~, n_frames] = size(rgbVid);
vid = rgbToGrayVid(rgbVid);

%% Read in Saliency Frames

saliencyMapHolder = zeros([18,18,448]);
%saliencyMapTYHolder = zeros([18,18,n_frames]);
for ii=1:5
    load(strcat('\ball_slowmo_crop_saliency_',num2str(ii),'.mat'));
    saliencyMapHolder(:,:,(ii-1)*100 + 1:min(ii*100,448)) = salMap(:,:,:);
    %saliencyMapTYHolder(:,:,(ii-1)*50 + 1:min(ii*50,n_frames)) = salMapTY(:,:,1:2:end);
    clear('regex','salMap*');
end
saliencyMapHolder = saliencyMapHolder(:,:,1:2:end-1);
size(saliencyMapHolder)

%% Compute optical flow frames
% Horn-Schunck Method
OF = opticalFlowHS();
flow = estimateFlow(OF, vid(:,:,1)); % Use 1st frame as reference
for i = 1:n_frames
    % I don't know how to initialize a vector of these objects
    flow(i) = estimateFlow(OF, vid(:,:,i));
end

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
    tmpSal=imresize(saliencyMapHolder(:,:,i),[rows, cols],'bilinear');
    imagesc(tmpSal);
    title('Saliency Map');
    
    subplot(133);
    imagesc(flow(i).Magnitude);
    title('Optical Flow');
    dur_calc = toc(st);
    pause(1/15 - dur_calc);
end
close(fig);
%% Compute OF Energy
of_energy = zeros(1,n_frames);
for i = 1:n_frames
    % Sum of all flow magnitudes - 'Overall movement'
    of_energy(i) = sum(sum(flow(i).Magnitude));
end
mov_avg_window = 15;

of_energy = movmean(of_energy, mov_avg_window);
of_energy = of_energy-min(of_energy);
of_energy = of_energy/max(of_energy);
%% Compute Saliency Energy
sal_energy = saliency_EMD_Avg_Energy(saliencyMapHolder);
%% Compare two graphs
of_fr = 4*(mean(of_energy)-of_energy);
of_fr = 60*2.^(of_fr);
sal_fr = 4*(mean(sal_energy)-sal_energy);
sal_fr = 60*2.^(sal_fr);

figure; 
subplot(121); hold on
plot(of_energy);
plot(sal_energy);
legend('Optical Flow','Saliency');
title('Energy Functions');

subplot(122); hold on
plot(of_fr);
plot(sal_fr);
legend('Optical Flow','Saliency');
title('Frame Rates');




%% Play Videos
figure, colormap('Gray');
currTime = 0;
for i = 1:n_frames
    imagesc(vid(:,:,i))
    title(['Original - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/60 );
    currTime = currTime + 1/60;
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





























