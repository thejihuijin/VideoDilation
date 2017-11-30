% RESIZE_VID Saves a video with new dimensions [newrows newcols]
% Input videos can be RGB or Greyscale
%
% inputName : Path to input video file
% outputName : Path and name for output video file
% newrows : Vertical size of resized video
% newcols : Horizontal size of resized video
function resize_vid(inputName,outputName,newrows,newcols)

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Read in Video
inputReader = VideoReader(inputName);

% Prepare output video
outputWriter = VideoWriter(outputName,'MPEG-4');
open(outputWriter);

while hasFrame(inputReader)
    % Read in Frame
	frame = readFrame(inputReader);
    % Resize frame and write
	outputFrame = imresize(frame, [newrows, newcols]);
	writeVideo(outputWriter, outputFrame);
	
end
close(outputWriter);
fprintf('done\n');
end

