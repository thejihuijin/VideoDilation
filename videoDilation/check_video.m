% CHECK_VIDEO Checks the dimension of the input video. The input video must
% be divisibl by 'dim', or the saliency algorithm will throw an error.
% If the dimensions don't fit the criteria, a resized video is generated
%
% input_video_path : path to input video
% dim : Dimension of video must be divisible by dim
%
% output_video_path : path to output video meeting size criteria

function [output_video_path] = check_video(input_video_path,dim)

% ECE6258: Digital image processing 
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
% Date Modified : 11/28/17
% By Erik Jorgensen (ejorgensen7@gatech.edu), Jihui Jin (jihui@gatech.edu)

% Check if video exists
if ~exist(input_video_path,'file')
    fprintf('Video does not exist');
    return;
end
% Read in video
inputReader = VideoReader(input_video_path);
rows = inputReader.Height;
cols = inputReader.Width;

% Check size
if mod(rows, dim) == 0 && mod(cols, dim) == 0
    output_video_path = input_video_path;
    return;
end

% Check if resized video exists
[filepath, name, ext] = fileparts(input_video_path);
output_video_path = strcat(filepath,'/rs_',name,ext);

if ~exist(output_video_path,'file') 
    newrows = dim*round(rows/dim);
    newcols = dim*round(cols/dim);
    fprintf('Resizing video to [%i %i]\n',newrows,newcols);
    resize_vid(input_video_path, output_video_path, newrows, newcols);
end
end

