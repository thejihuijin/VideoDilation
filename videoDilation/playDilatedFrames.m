function [] = playDilatedFrames( vidMat, frameIndices, fr, dilated_fr )
% PLAYDILATEDFRAMES Plays the frames designated by frameIndices within
% vidMat at the constant framerate fr

if length(size(vidMat)) == 4
    [~, ~, ~, frames] = size(vidMat);
else
    [~, ~, frames] = size(vidMat);
end

if frames < frameIndices(end)
    disp(['Frame index (' num2str(frameIndices(end)) ') greater than '...
        ' total number of frames (' num2str(frames) ')' ])
    frameIndices(end) = frames;
    disp(['Changed last frame index to ' num2str(frames)])
end

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