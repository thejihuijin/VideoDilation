# VideoDilation

This repo contains preliminary code for video dilation.
Video Dilation is an automated frame rate adjustment of videos 
depending on the energy present in each frame. More "interesting" events
are presented in slow motion while as less "interesting" events are
sped through.

## Contents
- Video Dilation : Library for automatically determining energy and adjusting frame rates
- Optical Flow : Optical Flow. Currently using Matlab Horn-Schunck
- Saliency : 3D FFT based saliency calculation*

## Dependencies
- Videos: Recommend short clips ~30 seconds long or less with minimal camera movement.
- Ran with Matlab_R2017b
  - Computer Vision Toolbox
  - Image Processing Toolbox
  
Example videos can be found here[link to be added later].
  
## How To Run
An example script is provided in `dilate_video.m`
Change the path to the video of interest at the top and run with the command:
```
dilate_video
```
in the matlab prompt.



























\* Saliency code provided by: Z. Long and G. AlRegib, "Saliency Detection for Videos Using 3D FFT Local Spectra," presented at Human Vision and Electronic Imaging XX, SPIE Electronic Imaging, San Francisco, CA, Feb. 8-12, 2015.
