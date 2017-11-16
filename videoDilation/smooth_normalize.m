function [smoothed_energy, std_dev] = smooth_normalize(energy, mov_avg_window)
%SMOOTH_NORMALIZE 
if ~exist('mov_avg_window')
    mov_avg_window=15;
end
std_dev = std(energy);

smoothed_energy=movmean(energy,mov_avg_window);
smoothed_energy = smoothed_energy-min(smoothed_energy);
smoothed_energy = smoothed_energy/max(smoothed_energy);

end

