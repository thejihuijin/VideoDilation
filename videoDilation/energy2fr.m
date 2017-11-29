function [fr_q_adj_smooth] = energy2fr(energy,fr,time_pad,scale)
%ENERGY2FR 

if ~exist('time_pad','var')
    time_pad=.2;
end
if ~exist('scale','var')
    scale=1;
end

energy_normal = 1-energy;

% Scale framerate to speed up
% Adjust framerate to exponential
fr_scaled = fr.*2.^(scale*(-1+2*energy_normal));


% time pad frame rate
fr_q_adj = adjustFR( fr_scaled, time_pad, fr );

% Smooth the adjusted framerate
mov_avg_window = 5;
fr_q_adj_smooth = movmean( fr_q_adj, mov_avg_window );
fr_q_adj_smooth = movmean( fr_q_adj_smooth, mov_avg_window );

end

