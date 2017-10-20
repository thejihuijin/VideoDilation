function grayVidMatrix = rgbToGrayVid( rgbVidMatrix )
%% Convert a 4D RGB matrix to a 3D grayscale matrix
%
% rgbVidMatrix : matrix (rows x cols x 3 x frames)
[rows, cols, ~, frames] = size(rgbVidMatrix);
grayVidMatrix = zeros(rows,cols,frames);

for i = 1:frames
    grayVidMatrix(:,:,i) = rgb2gray(rgbVidMatrix(:,:,:,i));
end

end

