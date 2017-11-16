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
[vid, fr] = vidToMat('gamecube04.mpg',500);
vid = rgbToGrayVid(vid);
[rows, cols, n_frames] = size(vid);
%% Read in Saliency Frames

saliencyMapHolder = zeros([rows/40,cols/40,n_frames]);
%saliencyMapTYHolder = zeros([18,18,n_frames]);
for ii=1:5
    load(strcat('data/saliency_mats/out\gamecube04_saliency_',num2str(ii),'.mat'));
    saliencyMapHolder(:,:,(ii-1)*100 + 1:min(ii*100,500)) = salMap(:,:,:);
    %saliencyMapTYHolder(:,:,(ii-1)*50 + 1:min(ii*50,n_frames)) = salMapTY(:,:,1:2:end);
    clear('regex','salMap*');
end
size(saliencyMapHolder)

%% Compute optical flow frames
% Horn-Schunck Method
OF = opticalFlowHS();
flow = estimateFlow(OF, vid(:,:,1)); % Use 1st frame as reference
for i = 1:n_frames
    % I don't know how to initialize a vector of these objects
    flow(i) = estimateFlow(OF, vid(:,:,i));
end
clear('OF');
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
    tmpSal=imresize(saliencyMapHolder(:,:,i),[rows, cols],'bilinear');
    imagesc(tmpSal,salminmax);
    title('Saliency Map');
    
    subplot(133);
    imagesc(flow(i).Magnitude,ofminmax);
    title('Optical Flow');
    dur_calc = toc(st);
    pause(1/15 - dur_calc);
end
close(fig);



%% Compute Standard Deviation
mario_stds = zeros(1,n_frames);
for ii=1:n_frames
    
    tmp = saliencyMapHolder(:,:,ii);
    mario_stds(ii) = std(tmp(:));
end
clear('tmp');
figure;
plot(mario_stds);
title('Standard Devations');




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

%% Compute masked Optical Flow energy

% Apply std threshold of .06
std_thresh = .06;
mask_energy_filtered = zeros(1,n_frames);

% Compute mask anyways
mask_energy_unfiltered = zeros(1,n_frames);

for ii=1:n_frames
    sal = imresize(saliencyMapHolder(:,:,ii),[rows, cols],'bilinear');
    sal_mask = imbinarize(sal);
    
    energy = sum(sum(sal_mask.*flow(ii).Magnitude));
    mask_energy_unfiltered(ii) = energy;
    
    tmp = saliencyMapHolder(:,:,ii);
    if std(tmp(:)) > std_thresh
        mask_energy_filtered(ii) = energy;
    end
    
end
clear('energy','tmp');

mask_energy_filtered = movmean(mask_energy_filtered,mov_avg_window);
mask_energy_filtered = mask_energy_filtered-min(mask_energy_filtered);
mask_energy_filtered = mask_energy_filtered/max(mask_energy_filtered);

mask_energy_unfiltered = movmean(mask_energy_unfiltered,mov_avg_window);
mask_energy_unfiltered = mask_energy_unfiltered-min(mask_energy_unfiltered);
mask_energy_unfiltered = mask_energy_unfiltered/max(mask_energy_unfiltered);
% figure; hold on
% plot(mask_energy_filtered)
% plot(mask_energy_unfiltered)
% legend('Filtered','Unfiltered')


%% Compare graphs
of_fr = 4*(mean(of_energy)-of_energy);
of_fr = 60*2.^(of_fr);
sal_fr = 4*(mean(sal_energy)-sal_energy);
sal_fr = 60*2.^(sal_fr);

figure; 
%subplot(121); 
hold on
plot(of_energy);
plot(sal_energy);
plot(mask_energy_filtered);
plot(mask_energy_unfiltered);
legend('Optical Flow','Saliency','Filtered Mask','Unfiltered Mask');
title('Energy Functions');

%subplot(122); hold on
%plot(of_fr);
%plot(sal_fr);
%legend('Optical Flow','Saliency');
%title('Frame Rates');




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





%% masked
unfiltered_fr = 4*(mean(mask_energy_unfiltered)-mask_energy_unfiltered);
unfiltered_fr = 60*2.^(unfiltered_fr);

figure, colormap('Gray')
currTime = 0;
for i = 1:n_frames
    imagesc(vid(:,:,i))
    title(['Saliency @ ',sprintf('%i',round(unfiltered_fr(i))),' fps - ' sprintf('%.2f',currTime) ' seconds elapsed'])
    pause( 1/unfiltered_fr(i) );
    currTime = currTime + 1/unfiltered_fr(i);
end

%% Compare Mask effect


frame = 250;
sal_img = imresize(saliencyMapHolder(:,:,frame),[rows, cols],'bilinear');
of_img = flow(frame).Magnitude;

figure;
subplot(221);
imagesc(sal_img)
title('Saliency Map');

subplot(222);
sal_mask = imbinarize(sal_img);
imagesc(sal_mask);
title('Binary Mask');

subplot(223);
imagesc(of_img);
title('Optical Flow');

subplot(224);
imagesc(of_img.*sal_mask);
title('Masked Optical Flow');

figure;
of_mask = imbinarize(of_img.*sal_mask);
imagesc(of_mask);



%%
of_sum = sum(sum(of_img.*sal_mask))
of_avg = mean(nonzeros(of_mask.*of_img))

mask_avg = zeros(1,n_frames);

for ii=1:n_frames
    sal = imresize(saliencyMapHolder(:,:,ii),[rows, cols],'bilinear');
    sal_mask = imbinarize(sal);
    
    of_img = flow(ii).Magnitude;
    of_mask = imbinarize(of_img.*sal_mask);
    
    energy = mean(nonzeros(of_mask.*of_img));
    mask_avg(ii) = energy;
end

mask_avg = mask_avg-min(mask_avg);
mask_avg = mask_avg/max(mask_avg);
mask_avg = movmean(mask_avg,mov_avg_window);
%% Compare graphs
of_fr = 4*(mean(of_energy)-of_energy);
of_fr = 60*2.^(of_fr);
sal_fr = 4*(mean(sal_energy)-sal_energy);
sal_fr = 60*2.^(sal_fr);

figure; 
%subplot(121); 
hold on
plot(of_energy);
plot(sal_energy);
plot(mask_energy_filtered);
plot(mask_energy_unfiltered);
plot(mask_avg);
legend('Optical Flow','Saliency','Filtered Mask','Unfiltered Mask','Masked Avg OF');
title('Energy Functions');

%subplot(122); hold on
%plot(of_fr);
%plot(sal_fr);
%legend('Optical Flow','Saliency');
%title('Frame Rates');



%% Compare total saliency

total_sal = sum(sum(saliencyMapHolder,1),2);
total_sal = total_sal(:);
total_sal = movmean(total_sal,mov_avg_window);
total_sal = total_sal-min(total_sal);
total_sal = total_sal/max(total_sal);
%% EMD from uniform PDF
edges = .025:.05:1;

%avgsalmap = mean(saliencyFrames,3);
%avg_hist = hist(avgsalmap(:),edges)/numel(avgsalmap);
%avg_cdf = cumsum(avg_hist);
emd_uni_energy = zeros([1 n_frames]);
for ii=1:n_frames
    salIm = saliencyMapHolder(:,:,ii);
    N = hist(salIm(:),edges)/numel(salIm);
    
    X = cumsum(N);
    emd_uni_energy(ii) = sum(abs(X-edges));
end
emd_uni_energy = smooth_normalize(emd_uni_energy);

%% Masked Saliency
masked_saliencyMapHolder = zeros(size(saliencyMapHolder));
for ii=1:n_frames
    temp = saliencyMapHolder(:,:,ii);
    masked_saliencyMapHolder(:,:,ii) = imbinarize(temp).*temp;
end
clear('temp');
masked_energy = saliency_EMD_Avg_Energy(masked_saliencyMapHolder);
%%
figure; hold on
%plot(total_sal,'DisplayName','total sal');
plot(of_energy,'DisplayName','total of');
plot(sal_energy,'DisplayName','sal emd');
plot(mask_avg,'DisplayName','Masked AVG OF');
%plot(emd_uni_energy,'DisplayName','uni emd');
plot(masked_energy,'DisplayName','masked emd avg');
title('Comparison of Energy Functions');
xlabel('Frame #')
ylabel('Energy')
legend('show');


%%
edges = .025:.05:1;
figure;
plot(emd_uni_energy)


%%

for ii=1:10
    imwrite(vid(:,:,ii*50),['mario_',num2str(ii*50),'.png'])
end

