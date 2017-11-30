% SLICEVID Convert a video file to a 4D matrix
% Assume input video is RGB 
% Dimensions = (rows, cols, 3, frames)
%
% INPUTS
% filename : String filename
% startTime : Time in video to start, in seconds
% endTime : Time in video to end, in seconds
% ds : Downsampling factor
%
% OUTPUTS
% vidMatrix : 4D matrix of video
% fr : framerate at which the video was taken (trunacted to an integer)

function [vidMatrix, fr] = sliceVid( filename, startTime, endTime, ds )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

v = VideoReader(filename);

rows = v.Height;
cols = v.Width;
fr = v.FrameRate;
fr = floor(fr);

total_frames = floor(v.Duration) * fr;

% Convert start and end times to frame numbers
start_frame = floor(startTime*fr)+1;
end_frame = floor(endTime*fr);
if end_frame > total_frames
    end_frame = total_frames;
end

num_frames = end_frame - start_frame;

% Downsampled vid size
rows_ds = floor(rows/ds);
cols_ds = floor(cols/ds);

vidMatrix = zeros(rows_ds, cols_ds, 3, num_frames);

% Throw away frames before start frame
for i = 1:start_frame-1
    readFrame(v);
end

% Read num_frames number of frames from the video
for i = 1:num_frames
    frame = readFrame(v);
    vidMatrix(:,:,:,i) = im2double(imresize(frame, [rows_ds cols_ds]));
end

end