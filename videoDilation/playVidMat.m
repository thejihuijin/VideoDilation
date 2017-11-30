% PLAYVIDMAT Plays a sequence of RGB or Grayscale frames with framerate 
% determined by input framerate vector or scalar.
% 
% INPUTS
% vidMat : 3D or 4D video matrix
% frameRates : Vector of variable framerates at which to play each frame of
% vidMat, or a constant at which to play the entire video.
% 
% OUTPUTS
% None - plays video in last opened figure

function [] = playVidMat(vidMat, frameRates)

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

% If input framerate is a constant, convert to a vector with that values
if length(frameRates) == 1
    frameRates = frameRates * ones(1,frames);
end

% Catch framerate vector length error and play until last available frame
if frames ~= length(frameRates)
    disp('Framerate vector and video matrix not equal length.')
    if length(frameRates) < frames
        frames = length(frameRates);
    end
end

% Play video frames sequentially by repeating calls to imagesc
currTime = 0;
if length(size(vidMat)) == 4
    for i = 1:frames
        imagesc(vidMat(:,:,:,i))
        axis off
        title([sprintf('%.2f',currTime) ' seconds elapsed'])
        pause( 1/frameRates(i) );
        currTime = currTime + 1/frameRates(i);
    end
else
    for i = 1:frames
        colormap gray
        imagesc(vidMat(:,:,i))
        title([sprintf('%.2f',currTime) ' seconds elapsed'])
        pause( 1/frameRates(i) );
        currTime = currTime + 1/frameRates(i);
    end
end      

end

