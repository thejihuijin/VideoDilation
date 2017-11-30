% QUANTFR Quantize the input framerate vector to a designated number of
% different values in the same range as the input vector. Depending on the
% rounding, this could create an extra quantized level which can simply be
% fixed by saturating values at the min/max of the input vector.
%
% INPUTS
% frVect : Variable framerate vector
% q_levels : Number of levels to quantize input vector into
% 
% OUTPUT
% frVect_quantized : Quantized version of input vector

function frVect_quantized = quantFR( frVect, q_levels )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Quantize vector to <q_levels> number of values within range of frVect
fr_range = range(frVect);
frVect_quantized = round(frVect/(fr_range/q_levels)) * (fr_range/q_levels);

end

