function [energy] = weight_pool(frames,p,mov_avg_window)
%WEIGHT_POOL Summary of this function goes here
%   Detailed explanation goes here
if ~exist('mov_avg_window','var')
    mov_avg_window=15;
end
if ~exist('p','var')
    p = 2;
end

[~,~,n_frames] = size(frames);
energy = zeros(1,n_frames);
for i = 1:n_frames
    tmpframe = frames(:,:,i);
    tmpframe = nonzeros(tmpframe(:));
    if ~isempty(tmpframe)
        energy(i) = sum(tmpframe.^(1+p))/sum(tmpframe.^p);
    end
end
energy = smooth_normalize(energy,mov_avg_window);
end

