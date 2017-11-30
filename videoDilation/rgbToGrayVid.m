% RGBTOGRAYVID Convert a 4D RGB matrix to a 3D grayscale matrix
%
% INPUT
% rgbVidMatrix : matrix (rows x cols x 3 x frames)
%
% OUTPUT
% 3D matrix of a grayscale video
function grayVidMatrix = rgbToGrayVid( rgbVidMatrix )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Sequentially convert each RGB frame to grayscale
[rows, cols, ~, frames] = size(rgbVidMatrix);
grayVidMatrix = zeros(rows,cols,frames);
for i = 1:frames
    grayVidMatrix(:,:,i) = rgb2gray(rgbVidMatrix(:,:,:,i));
end

end

