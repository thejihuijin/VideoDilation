function [vidMatrix, frame_rate] = sliceVid( filename, startTime, endTime, ds )
%% Convert a video file to a 4D matrix
% Assume video is RGB 
% Dimensions = (rows, cols, 3, frames)
%
% filename : String filename
v = VideoReader(filename);

rows = v.Height;
cols = v.Width;
frame_rate = v.FrameRate;
frame_rate = round(frame_rate);

total_frames = floor(v.Duration) * frame_rate;

start_frame = floor(startTime*frame_rate)+1;
end_frame = floor(endTime*frame_rate);
if end_frame > total_frames
    end_frame = total_frames;
end

num_frames = end_frame - start_frame;

% Downsampled vid size
rows_ds = floor(rows/ds);
cols_ds = floor(cols/ds);

vidMatrix = zeros(rows_ds, cols_ds, 3, num_frames);

% Throw away frames before
for i = 1:start_frame-1
    readFrame(v);
end
for i = 1:num_frames
    frame = readFrame(v);
%     vidMatrix(:,:,:,i) = im2double(frame(1:ds:end,1:ds:end,:));
    vidMatrix(:,:,:,i) = im2double(imresize(frame, [rows_ds cols_ds]));
end

end