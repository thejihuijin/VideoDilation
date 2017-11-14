%% 
clear all; close all; clc;

%% Load video
% [rgbVid, fr] = vidToMat('gamecube04.mpg');
% [rows, cols, ~, n_frames] = size(rgbVid);
% vid = rgbToGrayVid(rgbVid);


%% Load Subset of Video
nFrames = 200;
v = VideoReader('gamecube04.mpg');

rows = v.Height;
cols = v.Width;
dur = v.duration;
fr = v.FrameRate;
dur*fr
vidFrames = zeros([rows,cols, 3, nFrames]);
for ii = 1:nFrames
    frame = im2double(readFrame(v));
    vidFrames(:,:,:,ii) = frame;
end
%% load saliency data
nFrames = 200;
saliencyMapHolder = zeros([12,16,nFrames]);
saliencyMapTYHolder = zeros([12,16,nFrames]);
for ii=1:2
    load(strcat('../data/saliency_mats/out\gamecube04_saliency_',num2str(ii),'.mat'));
    saliencyMapHolder(:,:,(ii-1)*100 + 1:min(ii*100,nFrames)) = salMap(:,:,:);
    saliencyMapTYHolder(:,:,(ii-1)*100 + 1:min(ii*100,nFrames)) = salMapTY(:,:,:);
    clear('regex','salMap*');
end
%%
fig = figure;
for i = 1:nFrames
    st = tic;
    
    subplot(131);
    imagesc(vidFrames(:,:,:,i))
    %title('Constant framerate')
    
    subplot(132);
    tmpSal=imresize(saliencyMapHolder(:,:,i),[rows, cols],'bilinear');
    imagesc(tmpSal);
    title('Saliency Map');
    
    subplot(133);
    tmpSalT = imresize(saliencyMapTYHolder(:,:,i),[rows, cols],'bilinear');
    imagesc(tmpSalT);
    title('Temporal Saliency');
    
    dur_calc = toc(st);
    pause(.1 - dur_calc);
end
close(fig);
%% Calculate Metrics
% total
% average
% max

totalSaliency = reshape(sum(sum(saliencyMapHolder,1),2),[1 nFrames]);
meanSaliency = reshape(mean(mean(saliencyMapHolder,1),2),[1 nFrames]);
maxSaliency = reshape(max(max(saliencyMapHolder,[],1),[],2),[1 nFrames]);


%% EMD Plot for saliency maps
emd_avg = saliency_EMD_Avg_Energy(saliencyMapHolder);
%% Compare different metrics
figure; hold on;
plot(totalSaliency/max(totalSaliency));
plot(meanSaliency);
plot(maxSaliency);
%plot(emd_salmap/max(emd_salmap));
plot(emd_avg/max(emd_avg));
title('Comparison of Different Energy Functions')
legend('Total Saliency','Mean Saliency','maxSaliency','EMD vs avg frame');
%% Plot of EMD from Average Frame

figure;
plot(emd_avg);
xlabel('Frame #')
ylabel('EMD')
title('Smoothed EMD from avg Frame');

%% convert Video to Grey
grayFrames = zeros(rows, cols, nFrames);
for ii=1:nFrames
    grayFrames(:,:,ii) = rgb2gray(vidFrames(:,:,:,ii));
end
    
%% Play original, scaled, and smoothed in grayscale
clc;
f = figure, colormap('Gray')
for i = 1:nFrames
    imagesc(grayFrames(:,:,i))
    title(['Original @ ' num2str(fr) ' fps - ' ...
        sprintf('%.2f',i/fr) ' seconds elapsed'])
    pause( 1/fr );
end
close(f)
%% Adjust energy graph to frame rate
emd_frame_times = 2*(mean(emd_avg) - emd_avg);
figure; 
subplot(121); hold on;
plot(emd_frame_times)
plot(ones(size(emd_frame_times))*mean(emd_frame_times));
plot(ones(size(emd_frame_times))*max(emd_frame_times));
plot(ones(size(emd_frame_times))*min(emd_frame_times));

subplot(122);
plot(30*2.^(emd_frame_times));

emd_frame_times = 30*2.^(emd_frame_times);
%%
figure, colormap('Gray')
currTime = 0;
for i = 1:nFrames
    imagesc(grayFrames(:,:,i))
    title(['Variable framerate - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/emd_frame_times(i) );
    currTime = currTime + 1/emd_frame_times(i);
end

%%
figure(3), colormap('Gray')
currTime = 0;
for i = 1:nFrames
    imagesc(vid(:,:,i))
    title(['Smoothed variable framerate - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/fr_smoothed(i) );
    currTime = currTime + 1/fr_smoothed(i);
end

%% 
close all