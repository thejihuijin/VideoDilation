% SMOOTH_NORMALIZE filter and normalize a 1D array
% Assume input array is 1D
%
% filename : String filename
% startTime : Time in video to start, in seconds
% endTime : Time in video to end, in seconds
% ds : Downsampling factor
%
% vidMatrix : 4D matrix of video
% fr : framerate at which the video was taken (trunacted to an integer)

function [smoothed_energy, std_dev] = smooth_normalize(energy, mov_avg_window,mov_med_window)

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

%SMOOTH_NORMALIZE 
if ~exist('mov_avg_window','var')
    mov_avg_window=15;
end
if ~exist('mov_med_window','var')
    mov_med_window=5;
end
std_dev = std(energy);

smoothed_energy=movmedian(energy,mov_med_window);
smoothed_energy=movmean(smoothed_energy,mov_avg_window);
smoothed_energy = smoothed_energy-min(smoothed_energy);
smoothed_energy = smoothed_energy/max(smoothed_energy);

end

