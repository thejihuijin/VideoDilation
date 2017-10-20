%% Test video read and variable framerate
clear; clc; close all

v = VideoReader('traffic.avi');
rows = v.Height;
cols = v.Width;
dur = v.duration;
fr = v.FrameRate;
n_frames = floor(dur * fr);

% Put frames into 4D matrix (Rows, Cols, RGB, frames)
frames = zeros(rows, cols, 3, n_frames);
energy = zeros(n_frames,1);
for i = 1:n_frames
    im = readFrame(v);
    frames(:,:,:,i) = im;
    % Compute 'energy' as average pixel value
    energy(i) = mean(mean(rgb2gray(im)));
end
en_normal = energy - min(energy);
en_normal = en_normal / max(en_normal);

% Framerate normalized to range between 10 and 50
fr_scaled = en_normal*40 + 10;

% Play at original framerate
v = VideoReader('xylophone.mp4');
t_start = tic;
figure(1)
set(gcf, 'Position', [100, 600, 900, 350])
for i = 1:n_frames
    st = tic;
    vidFrame = readFrame(v);
    
    subplot(1,2,1)
    image(vidFrame)
    title('Constant framerate')
    
    subplot(1,2,2)
    plot(1:n_frames, fr*ones(n_frames,1), 'b-', i, fr, 'ro')
    title('Framerate per frame')
    xlabel('Frame #'), ylabel('Framerate')
    
    dur_calc = toc(st);
    pause( 1/fr - dur_calc );
end
toc(t_start)

% Play at energy-normalized framerate
t_start = tic;
figure(2)
set(gcf, 'Position', [100, 100, 900, 350])
for i = 1:n_frames
    st = tic;
    vidFrame = uint8(frames(:,:,:,i));
    
    subplot(1,2,1)
    image(vidFrame)
    title('Variable framerate')
    
    subplot(1,2,2)
    plot(1:n_frames, fr*ones(n_frames,1), 'k:', ...
        1:n_frames, fr_scaled, ...
        i, fr_scaled(i), 'ro')
    title('Framerate per frame')
    xlabel('Frame #'), ylabel('Framerate')
    
    dur_calc = toc(st);
    pause( 1/fr_scaled(i) - dur_calc );
end
toc(t_start)