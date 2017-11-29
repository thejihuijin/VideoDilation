% ADJUSTFR Extends the slowmo sections of a framerate vector to
% preemptively slow before 'exciting' segements and smoothly speed back up
% afterward.
%
% frVect : vector w/ time-varying framerate
% timeShift : Time to shift slow/speedups by, in seconds
% fr : Framerate video was taken at
%
% frAdjusted : vector w/ time-varying framerate & time-padded slow/speedup

function frAdjusted = adjustFR( frVect, timeShift, fr )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Compute derivative of FR
slope = conv(frVect, [1 -1], 'valid');

% Shift slowmo changes earlier, and shift speedup changes later
shift = round(timeShift*fr); 
slope_adj = slope;
for i = 1+shift : length(slope)-shift
    % Shift slowmo earlier
    if slope(i) < 0
        slope_adj(i-shift) = slope_adj(i-shift) + slope(i);
    % Shift speedup later
    elseif slope(i) > 0
        slope_adj(i+shift) = slope_adj(i+shift) + slope(i);
    end
    % Remove original slow/speedup
    slope_adj(i) = slope_adj(i) - slope(i);
end

% Recover framerate from derivative and reset range to min-max of frVect
frAdjusted = cumsum([0 slope_adj]) + frVect(1);
frAdjusted = min(frAdjusted, max(frVect));
frAdjusted = max(frAdjusted, min(frVect));

end

