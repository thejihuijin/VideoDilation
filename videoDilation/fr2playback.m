function playbackFrames = fr2playback( frameRates, playback_fr )
% Input a vector of framerates and a playback framerate
% Outputs a vector of frame indices to play sequentially

dilated_delays = 1 ./ frameRates;
dilated_times = [0 cumsum(dilated_delays(1:end))];

plybk_times = 0 : 1/playback_fr : dilated_times(end);% + 1/playback_fr;
playbackFrames = zeros(1,length(plybk_times));

frame_ptr = 1;
for i = 1:length(plybk_times)
    min_dist = inf;
    
    while 1
        if frame_ptr + 1 > length(dilated_times)
            if frame_ptr > length(playbackFrames)
                frame_ptr = length(playbackFrames);
            end
            playbackFrames(i) = frame_ptr;
            break
        end
        
        curr_dist = abs( dilated_times(frame_ptr) - plybk_times(i) );
        next_dist = abs( dilated_times(frame_ptr + 1) - plybk_times(i) );
        
        if curr_dist < min_dist
            min_dist = curr_dist;
        end
        
        if next_dist > curr_dist
            playbackFrames(i) = frame_ptr;
            break
        end
        
        frame_ptr = frame_ptr + 1;
    end
end

end

