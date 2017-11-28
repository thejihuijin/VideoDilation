function frAdjusted = adjustFR( frVect, time_shift, fr )

% Compute slope of FR
slope = conv(frVect, [1 -1], 'valid');

% Shift slowmos earlier, and speedups later
shift = round(time_shift*fr); 
slope_adj = slope;
for i = 1+shift : length(slope)-shift
    if slope(i) < 0
        slope_adj(i-shift) = slope_adj(i-shift) + slope(i);
        slope_adj(i) = 0;
    elseif slope(i) > 0
        slope_adj(i+shift) = slope_adj(i+shift) + slope(i);
        slope_adj(i) = 0;
    end
end

% Recover framerate from slope and reset range to min-max of frVect
frAdjusted = cumsum([0 slope_adj]) + frVect(1);
frAdjusted = min(frAdjusted, max(frVect));
frAdjusted = max(frAdjusted, min(frVect));

end

