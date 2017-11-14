function [] = playVidMat(vidMat, frameRates)
% PLAYVIDMAT Plays a sequence of RGB or Gray frames with framerate 
% determined by input framerate vector or scalar.

if length(size(vidMat)) == 4
    [~, ~, ~, frames] = size(vidMat);
else
    [~, ~, frames] = size(vidMat);
end

if length(frameRates) == 1
    frameRates = frameRates * ones(1,frames);
end
if frames ~= length(frameRates)
    error('Framerate vector and video matrix not equal length.')
end

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

