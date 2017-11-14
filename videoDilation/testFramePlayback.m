%% Test framerate quantization & frame interp/skipping
% ASSUMES testVidRead_3 HAS BEEN RUN ALREADY
clc

dilated_delays = 1 ./ fr_smoothed;
dilated_times = [0 cumsum(dilated_delays(1:end-1))];

plybk_fr = 100;
plybk_delays = 1 ./ (plybk_fr * ones(1,length(fr_smoothed)));
plybk_times = 0 : 1/plybk_fr : dilated_times(end) + 1/plybk_fr;

figure
plot(plybk_times, ones(1,length(plybk_times)), 'o', ...
    dilated_times, 1.1*ones(1,224), 'x')
axis([0 inf 0 2])

output_frames = zeros(1,length(plybk_times));
frame_ptr = 1;
for i = 1:length(plybk_times)
    closest_frame = 0;
    min_dist = inf;
    
    while 1
        curr_dist = abs( dilated_times(frame_ptr) - plybk_times(i) );
        next_dist = abs( dilated_times(frame_ptr + 1) - plybk_times(i) );
        
        if curr_dist < min_dist
            min_dist = curr_dist;
        end
        
        if next_dist > curr_dist
            output_frames(i) = frame_ptr;
            break
        end
        frame_ptr = frame_ptr + 1;
        
        if frame_ptr + 1 > length(dilated_times)
            break
        end
    end
end