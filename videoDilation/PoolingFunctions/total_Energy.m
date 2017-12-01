function [total_energy] = total_Energy(frames, mov_avg_window)
%TOTALENERGY Sums up energy per frame
if ~exist('mov_avg_window','var')
    mov_avg_window=15;
end
[~,~,n_frames] = size(frames);
total_energy = zeros(1,n_frames);
for i = 1:n_frames
    % Sum of all flow magnitudes - 'Overall movement'
    total_energy(i) = sum(sum(frames(:,:,i)));
end

total_energy = smooth_normalize(total_energy,mov_avg_window);
end
