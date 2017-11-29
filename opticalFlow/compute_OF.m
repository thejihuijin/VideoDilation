function [flow_mags] = compute_OF(vid)
%COMPUTE_OF Summary of this function goes here
%   Detailed explanation goes here

[rows,cols,n_frames] = size(vid);
OF = opticalFlowHS();
flow_mags = zeros(rows,cols,n_frames);
for i = 1:n_frames
    % Store magnitude frames
    flow = estimateFlow(OF, vid(:,:,i));
    flow_mags(:,:,i) = flow.Magnitude;
end


end

