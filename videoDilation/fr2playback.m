% FR2PLAYBACK Takes a variable framerate vector and finds the frames to be
% played at a constant framerate that best simulate the variable framerate.
%
% INPUTS
% frameRates : Vector of variable framerate per frame
% playback_fr : Constant framerate at which frames will be played
% 
% OUTPUT
% playbackFrames : Vector of frame indices into the video matrix to be 
% played sequentially at a constant FR, simulating a variable FR

function playbackFrames = fr2playback( frameRates, playback_fr )

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Delays between consecutive frames at variable framerate
dilated_delays = 1 ./ frameRates;

% Cumulative time until each frame at variable framerate
dilated_times = [0 cumsum(dilated_delays(1:end))];

% Cumulative time until each frame at constant framerate
plybk_times = 0 : 1/playback_fr : dilated_times(end);

% Find each frame from the variable framerate delays vector that occurs
% closest to each frame in the constant framerate delay vector.
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

