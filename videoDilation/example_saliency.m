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

%% EMD Histogram Comparison Example Single Frame
% sum(abs(cdf(x)-cdf(y)))
salIm = saliencyMapHolder(:,:,1);
figure;
subplot(221)
imagesc(salIm);
title('Saliency Image');
subplot(222)
histogram(salIm)
title('Histogram');

edges = .025:.05:1;
N = hist(salIm(:),edges)/numel(salIm);
subplot(223)
bar(edges,N)
title('PDF of Saliency Frame');

subplot(224)
uniform =ones(size(edges))/numel(salIm);
bar(edges,uniform);
title('PDF of Uniform')

figure;
X = cumsum(N);
Y = cumsum(uniform);
subplot(121)
plot(edges,X);
title('CDF of Saliency Frame');
subplot(122)
plot(edges,Y);
title('CDF of Uniform PDF');
fprintf('EMD = %f\n',sum(abs(X - Y)));
%% EMD Plot for saliency maps
emd_salmap = zeros([1 nFrames]);
edges = .025:.05:1;
uniform =ones(size(edges))/numel(saliencyMapHolder(:,:,1));
Y = cumsum(uniform);

avgsalmap = mean(saliencyMapHolder,3);
avg_hist = hist(avgsalmap(:),edges)/numel(avgsalmap);
avg_cdf = cumsum(avg_hist);
emd_avg = zeros([1 nFrames]);
for ii=1:nFrames
    salIm = saliencyMapHolder(:,:,ii);
    N = hist(salIm(:),edges)/numel(salIm);
    
    X = cumsum(N);
    emd_salmap(ii) = sum(abs(X-Y));
    emd_avg(ii) = sum(abs(X-avg_cdf));
end
%% Compare different metrics
figure; hold on;
plot(totalSaliency/max(totalSaliency));
plot(meanSaliency);
plot(maxSaliency);
plot(emd_salmap/max(emd_salmap));
plot(emd_avg/max(emd_avg));
title('Comparison of Different Energy Functions')
legend('Total Saliency','Mean Saliency','maxSaliency','Earth Movers Distance','EMD vs avg frame');
%% Plot of EMD from Average Frame
figure;
subplot(121);
plot(emd_avg)
xlabel('Frame #')
ylabel('EMD')
title('EMD from avg Frame');

subplot(122);
mov_avg_window = 15;
plot(smooth(emd_avg,10));
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
%%
mov_avg_window = 15;
emd_frame_times = smooth(emd_avg,mov_avg_window);
emd_frame_times = 1 - emd_frame_times/max(emd_frame_times);
figure; hold on

emd_frame_times = 5+50*emd_frame_times/max(emd_frame_times);
plot(emd_frame_times)
min_fr = 1/15;
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