# VideoDilation

This repo contains preliminary code for video dilation.
Video Dilation is an automated frame rate adjustment of videos 
depending on the energy present in each frame. More "interesting" events
are presented in slow motion while as less "interesting" events are
sped through.

## Authors
Jihui Jin and Erik Jorgensen

## Contents
- Video Dilation : Library for automatically determining energy and adjusting frame rates
- Optical Flow : Optical Flow. Currently using Matlab Horn-Schunck
- Saliency : 3D FFT based saliency calculation*

## Dependencies
- Videos: Recommend short clips ~30 seconds long or less with minimal camera movement.
- Ran with Matlab_R2017b
  - Computer Vision Toolbox
  - Image Processing Toolbox
  
Example videos can be found [here](https://www.dropbox.com/sh/wpze1o1taqz6yyh/AABEjfNdWdxFotm40nC9Dp_ma?dl=0 "Test Videos").
  
## How To Run
An example script is provided in `dilate_video.m`
Change the path to the video of interest at the top and run with the command:
```
dilate_video
```
in the matlab prompt.


The last section includes code to save the resulting video to a file with the option of a filename.
```
% Retrieves original rgb frames
rgbvid = vidToMat(filename);

% saves as output 'test.mp4' in current directory
saveDilatedFrames( rgbvid , adjusted_smooth_playback, fr, time_padded_fr);

% Can set optional output directory/filename
saveDilatedFrames( rgbvid , adjusted_smooth_playback, fr, time_padded_fr, outputfilename )
```

### Note
There are minor potential issues with the video resize code. All videos included in the dropbox
work on Mac MATLAB_R2017b, but the resizing has caused some distortion on a windows machine. This 
problem does not affect every video. We have isolated the issue to the readframe call, but have not 
fixed the issues since then. 






















\* Saliency code provided by: Z. Long and G. AlRegib, "Saliency Detection for Videos Using 3D FFT Local Spectra," presented at Human Vision and Electronic Imaging XX, SPIE Electronic Imaging, San Francisco, CA, Feb. 8-12, 2015.
