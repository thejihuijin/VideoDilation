%% Test video read and variable framerate - #2
% 'rhinos.avi'
% 'xylophone.mp4'
% 'traffic.avi'
% 'visiontraffic.avi'
clear; clc; close all

%% Load video
[rgbVid, fr] = vidToMat('traffic.avi');
[rows, cols, ~, n_frames] = size(rgbVid);
vid = rgbToGrayVid(rgbVid);

%% Compute optical flow frames
% Lucas Kanade Method
% get coordinate for u and v in the original frame
[X,Y] = meshgrid(1:cols, 1:rows);
X_ds = X(1:10:end, 1:10:end);
Y_ds = Y(1:10:end, 1:10:end);

u_vid = zeros(size(vid));
v_vid = zeros(size(vid));
u_vid_ds = zeros(rows/10,cols/10,n_frames);
v_vid_ds = zeros(rows/10,cols/10,n_frames);
tstart = tic;
for i = 2:n_frames
    [u, v] = optFlowLK(vid(:,:,i-1), vid(:,:,i), 40);
    
    u_vid(:,:,i) = u; % Horizontal flow magnitudes
    v_vid(:,:,i) = v; % Vertical flow magnitudes

    % Downsample flow vectors
    u_vid_ds(:,:,i) = u_vid(1:10:end, 1:10:end, i);
    v_vid_ds(:,:,i) = v_vid(1:10:end, 1:10:end, i);
    
    if mod(i,5) == 0
        disp(['Compute time up to frame ' num2str(i) ...
            ': ' sprintf('%.2f',toc(tstart))])
    end
end

%% Compute energy vector
energy = zeros(n_frames,1);
flow_mag = zeros(size(u_vid));
for i = 1:n_frames
    flow_mag(:,:,i) = sqrt(u_vid(:,:,i).^2 + v_vid(:,:,i).^2);
    energy(i) = sum(sum(flow_mag(:,:,i)));
end
    
    
%% Normalize energy vector
energy_normal = (energy - min(energy)) / (max(energy - min(energy)));

% Invert energy scale so greater energy correlates to lower framerate
energy_normal = 1 - energy_normal;

%% Scale framerate to desired range
% Framerate normalized to range between 10 and 30
range_fr = 20;
min_fr = 10;
fr_scaled = energy_normal*range_fr + min_fr;

%% Smooth, renormalize, and scale the energy-normalized framerate
mov_avg_window = 15;
energy_smoothed = movmean(energy_normal, mov_avg_window);
energy_smoothed = (energy_smoothed - min(energy_smoothed)) ...
    / (max(energy_smoothed - min(energy_smoothed)));

fr_smoothed = energy_smoothed*range_fr + min_fr;

%% Play at original framerate
sf = 0.5;
slow_mo = 3;

figure(1)
set(gcf, 'Position', [100, 600, 900, 350])
for i = 1:n_frames
    st = tic;
    
    subplot(1,2,1)
    imagesc(rgbVid(:,:,:,i))
    title('Constant framerate')
    hold on
    quiver(X_ds, Y_ds, sf*u_vid_ds(:,:,i), sf*v_vid_ds(:,:,i), ...
        'color', [0 1 0], 'AutoScale', 'off')
    hold off
    
    subplot(1,2,2)
    plot(1:n_frames, fr*ones(n_frames,1), 'b-', i, fr, 'ro')
    title([sprintf('%.2f',i/fr) ' seconds elapsed'])
    xlabel('Frame #'), ylabel('Framerate')
    
    pause_dur = 1/fr - toc(st);
    pause( slow_mo*pause_dur );
end

%% Play at energy-normalized framerate
figure(2)
set(gcf, 'Position', [100, 350, 900, 350])
currTime = 0;
for i = 1:n_frames
    st = tic;
    
    subplot(1,2,1)
    imagesc(rgbVid(:,:,:,i))
    title('Variable framerate')
    hold on
    quiver(X_ds, Y_ds, sf*u_vid_ds(:,:,i), sf*v_vid_ds(:,:,i), ...
        'color', [0 1 0], 'AutoScale', 'off')
    hold off
    
    subplot(1,2,2)
    plot(1:n_frames, fr*ones(n_frames,1), 'k:', ...
        1:n_frames, fr_scaled, ...
        i, fr_scaled(i), 'ro')
    title([sprintf('%.2f',currTime) ' seconds elapsed'])
    xlabel('Frame #'), ylabel('Framerate')
    
    currTime = currTime + 1/fr_scaled(i);
    pause_dur = 1/fr_scaled(i) - toc(st);
    pause( slow_mo*pause_dur );
end

%% Play at smoothed-energy-normalized framerate
figure(3)
set(gcf, 'Position', [100, 100, 900, 350])
currTime = 0;
for i = 1:n_frames
    st = tic;
    
    subplot(1,2,1)
    imagesc(rgbVid(:,:,:,i))
    title('Variable framerate')
    hold on
    quiver(X_ds, Y_ds, sf*u_vid_ds(:,:,i), sf*v_vid_ds(:,:,i), ...
        'color', [0 1 0], 'AutoScale', 'off')
    hold off
    
    subplot(1,2,2)
    plot(1:n_frames, fr*ones(n_frames,1), 'k:', ...
        1:n_frames, fr_smoothed, ...
        i, fr_smoothed(i), 'ro')
    title([sprintf('%.2f',currTime) ' seconds elapsed'])
    xlabel('Frame #'), ylabel('Framerate')
    
    currTime = currTime + 1/fr_smoothed(i);
    pause_dur = 1/fr_smoothed(i) - toc(st);
    pause( slow_mo*pause_dur );
end