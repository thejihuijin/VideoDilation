function [energy] = five_num_sum(frames,mov_avg_window)
%FIVE_NUM_SUM Summary of this function goes here
%   Detailed explanation goes here
if ~exist('mov_avg_window','var')
    mov_avg_window=15;
end
[~,~,n_frames] = size(frames);
energy = zeros(1,n_frames);
for i = 1:n_frames
    tmpframe = frames(:,:,i);
    tmpframe = nonzeros(tmpframe(:));
    if ~isempty(tmpframe)
        energy(i) = mean(tmpframe) + max(tmpframe)...
                + sum(prctile(tmpframe,[25 50 75]));
        energy(i) = energy(i)/5;
    end
end

energy = smooth_normalize(energy,mov_avg_window);
end


