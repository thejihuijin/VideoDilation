% SAVEDILATEDFRAMES saves the frames as designated by the vector of
% indices, frameIndices, at a constant framerate.
% 
% INPUTS
% vidMat : 3D or 4D video matrix
% frameIndices : Vector of indices into vidMat to be played sequentially
% fr : Constant framerate at which to play frames
% dilated_fr : Variable framerate at which each frame from vidMat is played
% outputName : Path to output video
% 
% OUTPUTS
% None - saves dilated video 

function [] = playDilatedFrames( vidMat, frameIndices, fr, dilated_fr, outputName )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Grab number of frames
if length(size(vidMat)) == 4
    [~, ~, ~, frames] = size(vidMat);
else
    [~, ~, frames] = size(vidMat);
end

% Check outputName
if ~exist('outputName','var')
    outputName = 'test.mp4';
end

% Catch off-by-one error and correct
if frames < frameIndices(end)
    disp(['Frame index (' num2str(frameIndices(end)) ') greater than '...
        ' total number of frames (' num2str(frames) ')' ])
    frameIndices(end) = frames;
    disp(['Changed last frame index to ' num2str(frames)])
end

% Play video frames sequentially by repeating calls to imagesc
outputWriter = VideoWriter(outputName,'MPEG-4');
outputWriter.FrameRate = fr;
open(outputWriter);

if length(size(vidMat)) == 4
    for i = 1:length(frameIndices)
        writeVideo(outputWriter, vidMat(:,:,:,frameIndices(i)));
    end
else
    for i = 1:length(frameIndices)
        writeVideo(outputWriter, vidMat(:,:,frameIndices(i)));
    end
end      
fprintf('Video saved to %s\n',outputName);
end