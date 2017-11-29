function [outputArg1,outputArg2] = resize_vid(inputName,outputName,newrows,newcols)
inputReader = VideoReader(inputName);

outputWriter = VideoWriter(outputName,'MPEG-4');
open(outputWriter);

while hasFrame(inputReader)
	frame = readFrame(inputReader);
	outputFrame = imresize(frame, [newrows, newcols]);
	writeVideo(outputWriter, outputFrame);
	
end
close(outputWriter);
fprintf('done\n');
end

