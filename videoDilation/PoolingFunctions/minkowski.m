function [energy] = minkowski(frames, p, mov_avg_window)
%MINKOWSKI Summary of this function goes here
%   Detailed explanation goes here
if ~exist('mov_avg_window','var')
    mov_avg_window=15;
end
if ~exist('p','var')
    p = 2;
end
[~,~,n_frames] = size(frames);
frames = frames.^p;
energy = zeros(1,n_frames);
for i=1:n_frames
    tmp = nonzeros(frames(:,:,i));
    if ~isempty(tmp)
        energy(i) = mean(tmp(:));
    end
end
energy = smooth_normalize(energy,mov_avg_window);

end

