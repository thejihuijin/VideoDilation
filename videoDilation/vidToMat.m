function [vidMatrix, frame_rate] = vidToMat( filename, max_frames )
%% Convert a video file to a 4D matrix
% Assume video is RGB 
% Dimensions = (rows, cols, 3, frames)
%
% filename : String filename
% max_frames: int max frames to read in
v = VideoReader(filename);

rows = v.Height;
cols = v.Width;
frame_rate = v.FrameRate;
num_frames = floor(v.duration * frame_rate);
if exist('max_frames','var')
     % third parameter does not exist, so default it to something
      num_frames = min(max_frames,num_frames);
end
vidMatrix = zeros(rows, cols, 3, num_frames);
for i = 1:num_frames
    frame = readFrame(v);
    vidMatrix(:,:,:,i) = im2double(frame);
end

end

