% PLAYDILATEDFRAMES Plays the frames as designated by the vector of
% indices, frameIndices, at a constant framerate.
% 
% INPUTS
% vidMat : 3D or 4D video matrix
% frameIndices : Vector of indices into vidMat to be played sequentially
% fr : Constant framerate at which to play frames
% dilated_fr : Variable framerate at which each frame from vidMat is played
% 
% OUTPUTS
% None - plays video in last opened figure

function [] = playDilatedFrames( vidMat, frameIndices, fr, dilated_fr )

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

% Catch off-by-one error and correct
if frames < frameIndices(end)
    disp(['Frame index (' num2str(frameIndices(end)) ') greater than '...
        ' total number of frames (' num2str(frames) ')' ])
    frameIndices(end) = frames;
    disp(['Changed last frame index to ' num2str(frames)])
end

% Play video frames sequentially by repeating calls to imagesc
currTime = 0;
if length(size(vidMat)) == 4
    for i = 1:length(frameIndices)
        imagesc(vidMat(:,:,:,frameIndices(i)))
        axis off
        title([sprintf('%.2f',currTime) ' seconds elapsed at ' ...
            sprintf('%.1f',dilated_fr(frameIndices(i))) ' frames/s'])
        pause( 1/fr );
        currTime = currTime + 1/fr;
    end
else
    for i = 1:length(frameIndices)
        colormap gray
        imagesc(vidMat(:,:,frameIndices(i)))
        title([sprintf('%.2f',currTime) ' seconds elapsed at ' ...
            sprintf('%.1f', dilated_fr(frameIndices(i))) ' frames/s'])
        pause( 1/fr );
        currTime = currTime + 1/fr;
    end
end      

end