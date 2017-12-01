
% VIDTOMAT Convert a RGB video file to a 4D matrix of image frames
%
% INPUT
% filename : String filename of video
% 
% OUTPUT
% vidMatrix : 4D matrix of rgb frames
% frame_rate : Framerate at which the video file was encoded

function [vidMatrix, frame_rate] = vidToMat( filename )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)


v = VideoReader(filename);

rows = v.Height;
cols = v.Width;
frame_rate = v.FrameRate;
num_frames = floor(v.duration * frame_rate);

% Read each frame sequentially
vidMatrix = zeros(rows, cols, 3, num_frames);
for i = 1:num_frames
    frame = readFrame(v);
    vidMatrix(:,:,:,i) = im2double(frame);
end

end

