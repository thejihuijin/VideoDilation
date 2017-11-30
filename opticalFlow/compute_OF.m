% COMPUTE_OF returns the optical flow magnitudes for a given video
% Currently uses Horn Schunck
% Assumes input video is in grey scale
% Dimensions = (rows, cols, frames)
%
% vid : 3D video matrix
%
% flow_mags : 3D Optical Flow magnitudes

function [flow_mags] = compute_OF(vid)

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

[rows,cols,n_frames] = size(vid);
OF = opticalFlowHS();
flow_mags = zeros(rows,cols,n_frames);
for i = 1:n_frames
    % Store magnitude frames
    flow = estimateFlow(OF, vid(:,:,i));
    flow_mags(:,:,i) = flow.Magnitude;
end


end

