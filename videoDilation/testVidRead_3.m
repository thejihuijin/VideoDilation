%% Test video read and variable framerate - #2
clear; clc; close all
% 'rhinos.avi'
% 'xylophone.mp4'
% 'traffic.avi'
% 'visiontraffic.avi'
% [rgbVid, fr] = vidToMat('traffic.avi');
% [rgbVid, fr] = sliceVid('go_pro.mp4', 17, 28, 5);
% [rgbVid, fr] = sliceVid('rockBeach.mp4', 18, 23, 5);
% [rgbVid, fr] = sliceVid('ball_slowmo.mp4', 0, 13, 10);
% [rgbVid, fr] = sliceVid('ball_slowmo_crop.mp4', 0, 15, 2);
% [rgbVid, fr] = sliceVid('../../ball_slowmo_crop.mp4', 0, 15, dim_ds);
% [rgbVid, fr] = sliceVid('../../S9CEzbd7Ftc.mp4', 30, 40, dim_ds);

%% Load video
dim_ds = 4; % Downsampling of height & width dimensions
[rgbVid, fr] = sliceVid('../../ball_slowmo_crop.mp4', 0, 15, dim_ds);

% rgbVid = cat(4, rgbVid, flip(rgbVid,4), rgbVid);

slow_speed = 1/8; % If video is in slow-mo change this var
fr = fr/slow_speed;

frame_ds = 1; % Downsample the number of frames
rgbVid = rgbVid(:,:,:,1:frame_ds:end);
fr = fr/frame_ds;

% Convert to gray & record final dimensions
[rows, cols, ~, n_frames] = size(rgbVid);
vid = rgbToGrayVid(rgbVid);

clear rgbVid

%% Compute optical flow frames
% Horn-Schunck Method
OF = opticalFlowHS();
flow = estimateFlow(OF, vid(:,:,1)); % Set 1st frame as reference
for i = 1:n_frames
    % I don't know how to initialize a vector of these objects
    flow(i) = estimateFlow(OF, vid(:,:,i));
end

%% Compute energy vector
energy = zeros(1,n_frames);
for i = 1:n_frames
    % Sum of all flow magnitudes - 'Overall movement'
    energy(i) = sum(sum(flow(i).Magnitude));
end

% Interpolate outliers with median
energy = movmedian(energy,5);

%% Normalize energy vector
energy_normal = (energy - min(energy)) / range(energy);

% Invert energy scale so greater energy correlates to lower framerate
energy_normal = 1 - energy_normal;

%% Scale framerate to desired range
% Framerate normalized to range between min_fr and min_fr + range_fr
min_fr = 120;
range_fr = 240;
fr_scaled = energy_normal*range_fr + min_fr;

% Alternatively, scale so that the mean framerate = original framerate
energy_normal2 = energy_normal / mean(energy_normal);
fr_scaled2 = min_fr + (fr - min_fr)*energy_normal2;

%% Smooth, renormalize, and scale the energy-normalized framerate
mov_avg_window = 15;

energy_smooth = movmean(energy_normal, mov_avg_window);
energy_smooth = (energy_smooth - min(energy_smooth))/range(energy_smooth);
fr_smoothed = energy_smooth*range_fr + min_fr;
% fr_smoothed = movmean(fr_scaled, mov_avg_window);

% Alternatively, scale so that the mean framerate = original framerate
energy_smoothed2 = movmean( energy_normal, mov_avg_window);
energy_smoothed2 = energy_smoothed2 / mean(energy_smoothed2);
fr_smoothed2 = min_fr + (fr - min_fr)*energy_smoothed2;

%% Quantize framerate, adjust framerate changes, re-smooth framerate
q_levels = 5;
fr_quant = quantFR( fr_smoothed, q_levels );

fr_q_adj = adjustFR( fr_quant, 0.2, fr );

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

%% Play original, scaled, and smoothed in grayscal
figure
playVidMat( vid, fr )
% 
% figure(2)
% playVidMat( vid, fr_scaled )
% 
% figure(3)
% playVidMat( vid, fr_smoothed )

figure
playDilatedFrames( vid, scaled_playback, fr, fr_scaled )

figure
playDilatedFrames( vid, smoothed_playback, fr, fr_smoothed )

figure
playDilatedFrames( vid, adjusted_smooth_playback, fr, fr_q_adj_smooth )

% %% Play at original framerate
% slow_mo = 1; % Number of times slower to play videos
% 
% figure(1)
% set(gcf, 'Position', [100, 600, 900, 350])
% for i = 1:n_frames
%     st = tic;
%     
%     % Plot original video with optical flow vectors
%     subplot(1,2,1)
%     imagesc(rgbVid(:,:,:,i))
%     title('Constant framerate')
%     hold on
%     plot(flow(i), 'DecimationFactor', [5 5], 'ScaleFactor', 30)
%     hold off
%     
%     subplot(1,2,2)
%     plot(1:n_frames, fr*ones(n_frames,1), 'b-', i, fr, 'ro')
%     title([sprintf('%.2f',i/fr) ' seconds elapsed'])
%     xlabel('Frame #'), ylabel('Framerate')
%     
%     pause_dur = 1/fr - toc(st);
%     pause( slow_mo*pause_dur );
% end
% 
% %% Play at energy-normalized framerate
% figure(2)
% set(gcf, 'Position', [100, 350, 900, 350])
% currTime = 0;
% for i = 1:n_frames
%     st = tic;
%     
%     % Plot energy-normed video with optical flow vectors
%     subplot(1,2,1)
%     imagesc(rgbVid(:,:,:,i))
%     title('Variable framerate')
%     hold on
%     plot(flow(i), 'DecimationFactor', [5 5], 'ScaleFactor', 30)
%     hold off
%     
%     subplot(1,2,2)
%     plot(1:n_frames, fr*ones(n_frames,1), 'k:', ...
%         1:n_frames, fr_scaled, ...
%         i, fr_scaled(i), 'ro')
%     title([sprintf('%.2f',currTime) ' seconds elapsed'])
%     xlabel('Frame #'), ylabel('Framerate')
%     
%     currTime = currTime + 1/fr_scaled(i);
%     pause_dur = 1/fr_scaled(i) - toc(st);
%     pause( slow_mo*pause_dur );
% end
% 
% %% Play at smoothed-energy-normalized framerate
% figure(3)
% set(gcf, 'Position', [100, 100, 900, 350])
% currTime = 0;
% for i = 1:n_frames
%     st = tic;
%     
%     % Plot smoothed, energy-normed video with optical flow vectors
%     subplot(1,2,1)
%     imagesc(rgbVid(:,:,:,i))
%     title('Variable framerate')
%     hold on
%     plot(flow(i), 'DecimationFactor', [5 5], 'ScaleFactor', 30)
%     hold off
%     
%     subplot(1,2,2)
%     plot(1:n_frames, fr*ones(n_frames,1), 'k:', ...
%         1:n_frames, fr_smoothed, ...
%         i, fr_smoothed(i), 'ro')
%     title([sprintf('%.2f',currTime) ' seconds elapsed'])
%     xlabel('Frame #'), ylabel('Framerate')
%     
%     currTime = currTime + 1/fr_smoothed(i);
%     pause_dur = 1/fr_smoothed(i) - toc(st);
%     pause( slow_mo*pause_dur );
% end