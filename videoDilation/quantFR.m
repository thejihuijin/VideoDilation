function frVect_quantized = quantFR( frVect, q_levels )

fr_range = range(frVect);

% Quantize the input vector to q_levels different values in the same range
frVect_quantized = round(frVect/(fr_range/q_levels)) * (fr_range/q_levels);

end

