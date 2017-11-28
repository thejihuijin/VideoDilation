function frVect_quantized = quantFR( frVect, q_levels )
% QUANTFR Quantize the input framerate vector to <q_levels> different 
% values in the same range as the input vector
fr_range = range(frVect);

frVect_quantized = round(frVect/(fr_range/q_levels)) * (fr_range/q_levels);

end

