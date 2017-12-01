% ENERGY2FR Converts energy to time padded frame rate
% Inverts an energy function, then scales it between -scale to scale and
% passes it through an exponential to determine speed up factor.
% Frame rate is then padded to begin slow down prior to "interesting"
% events
%
% energy : 1D energy function. High values will be slowed down
% fr : Original frame rate 
% time_pad : Time to shift slows/speedups by, in seconds
% scale : Speed up/slow down factor, 2^scale
%
% fr_adj_smooth : new frame rate for each frame in original video

function [fr_adj_smooth] = energy2fr(energy,fr,time_pad,scale)

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

if ~exist('time_pad','var')
    time_pad=.2;
end
if ~exist('scale','var')
    scale=1;
end

energy_normal = 1-energy;

% Scale framerate to speed up
% Adjust framerate to exponential
fr_scaled = fr.*2.^(scale*(2*(energy_normal-mean(energy_normal))));


% time pad frame rate
fr_adj = adjustFR( fr_scaled, time_pad, fr );

% Smooth the adjusted framerate
mov_avg_window = 5;
fr_adj_smooth = movmean( fr_adj, mov_avg_window );
fr_adj_smooth = movmean( fr_adj_smooth, mov_avg_window );

end

