function [saliencyMap, saliencyMapTime] = compute_saliency(iFName,wsize_s,wsize_t,...
    wsize_fs,wsize_ft,scaleFactor,segLength)
%
% Generate saliency map and time weighted saliency map for given video
%
% iFName: input video file name
% wsize_s: spatial window size (x and y) for center-surround comparison
% wsize_t: temporal window size (t) for c-s comparison
% wsize_fs: spatial window size (x and y) for fft
% wsize_ft: temporal window size (t) for fft
% scaleFactor: factor for scaling down video frames spatially
% segLength: # of frames to be processed at one time
%
% saliencyMap: saliency maps for every frame in the video
% saliencyMapTime: time weighted saliency maps for every frame in the video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Please cite our paper when you use the code:
%
% Z. Long and G. AlRegib, "Saliency Detection for Videos Using 3D FFT Local Spectra," presented at Human Vision and Electronic Imaging XX, SPIE Electronic Imaging, San Francisco, CA, Feb. 8-12, 2015.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCLAIMER: The Matlab code provided is only for evaluation of the algorithm. 
% Neither the authors of the code, nor affiliations of the authors can be held 
% responsible for any damages arising out of using this code in any manner. 
% Please use the code at your own risk.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of program
% Parameters
hw_t=(wsize_t-1)/2; % Half window size (t) for c-s comparison
hw_ft=(wsize_ft-1)/2; % Half window size (t) for fft

% Prepare input/output video files
fp=VideoReader(iFName);
frameH=fp.Height; frameW=fp.Width;
frameHS=frameH*scaleFactor; frameWS=frameW*scaleFactor;

% Read and process the frames
count0=0; count=0; nSeg=0;
jStop=fp.NumberOfFrames;

% Generate holders for output
saliencyMap = zeros(frameHS/wsize_fs,frameWS/wsize_fs,jStop);
saliencyMapTime = zeros(frameHS/wsize_fs,frameWS/wsize_fs,jStop);
index = 0;
tweight = 3/4;
for j=1:jStop
    frame=read(fp,j);
    frameYIQ=rgb2ntsc(frame);
    count=count+1;

    frameY(:,:,count)=single(frameYIQ(:,:,1)); % For results display
    frameI(:,:,count)=single(frameYIQ(:,:,2)); % For results display
    frameQ(:,:,count)=single(frameYIQ(:,:,3)); % For results display

    tmp=imresize(frameYIQ(:,:,1),[frameHS frameWS],'bilinear');
    stdTmp=std(tmp(:));
    if stdTmp==0 % Normalize, make sure not to divide by zero
        display(['detSal.m: Y: j=' num2str(j) ', std=0']);
        dataY(:,:,count)=single(tmp-mean(tmp(:)));
    else
        dataY(:,:,count)=single((tmp-mean(tmp(:)))/stdTmp);
    end
    
    tmp=imresize(frameYIQ(:,:,2),[frameHS frameWS],'bilinear');
    stdTmp=std(tmp(:));
    if stdTmp==0 % Normalize, make sure not to divide by zero
        display(['detSal.m: I: j=' num2str(j) ', std=0']);
        dataI(:,:,count)=single(tmp-mean(tmp(:)));
    else
        dataI(:,:,count)=single((tmp-mean(tmp(:)))/stdTmp);
    end
    
    tmp=imresize(frameYIQ(:,:,3),[frameHS frameWS],'bilinear');
    stdTmp=std(tmp(:));
    if stdTmp==0 % Normalize, make sure not to divide by zero
        display(['detSal.m: Q: j=' num2str(j) ', std=0']);
        dataQ(:,:,count)=single(tmp-mean(tmp(:)));
    else
        dataQ(:,:,count)=single((tmp-mean(tmp(:)))/stdTmp);
    end
    
    if (count-count0)==(segLength+hw_t+hw_ft)... % Extra data at border for accuracy
            || j==jStop
        % Calculate FFT
        tic;
        dataFY=applyFft(dataY,wsize_fs,wsize_ft);
        dataFI=applyFft(dataI,wsize_fs,wsize_ft);
        dataFQ=applyFft(dataQ,wsize_fs,wsize_ft);
        disp(['Calulating FFT: ' num2str(toc) 'sec']);
        
        % Compute space-time saliency
        tic;
        wS=0.5; wT=0.5; % To be improved with adaptive weights
        [salMapY,salMapTY,salMapSY]=calcSalMap(dataFY,wsize_s,wsize_t,wS,wT); clear dataFY;
        [salMapI,salMapTI,salMapSI]=calcSalMap(dataFI,wsize_s,wsize_t,wS,wT); clear dataFI;
        [salMapQ,salMapTQ,salMapSQ]=calcSalMap(dataFQ,wsize_s,wsize_t,wS,wT); clear dataFQ;
        salMap=(salMapY+salMapI+salMapQ)/3.0; % To be improved with adaptive weights
        disp(['Calculating saliency map: ' num2str(toc) 'sec']);

        % Save saliency map for the segment
        if (j==jStop)
            saliencyMap(:,:,index*segLength+1:jStop) = salMap(:,:,count0+1:count);
            
            salMapTime = 1/4*((salMapSI + salMapSQ + salMapSY)/3) ...
                       + 3/4*((salMapTI + salMapTQ + salMapTY)/3);

            saliencyMapTime(:,:,index*segLength+1:jStop) = salMapTime(:,:,count0+1:count);

%            salMap=salMap(:,:,count0+1:count);
%             salMapY=salMapY(:,:,count0+1:count);
%             salMapI=salMapI(:,:,count0+1:count);
%             salMapQ=salMapQ(:,:,count0+1:count);
%             salMapTY=salMapTY(:,:,count0+1:count);
%             salMapTI=salMapTI(:,:,count0+1:count);
%             salMapTQ=salMapTQ(:,:,count0+1:count);
%             salMapSY=salMapSY(:,:,count0+1:count);
%             salMapSI=salMapSI(:,:,count0+1:count);
%             salMapSQ=salMapSQ(:,:,count0+1:count);
        else
            saliencyMap(:,:,index*segLength+1:(index+1)*segLength) = salMap(:,:,count0+1:count0+segLength);
            salMapTime = (1-tweight)*((salMapSI + salMapSQ + salMapSY)/3) ...
                       + tweight*((salMapTI + salMapTQ + salMapTY)/3);
            saliencyMapTime(:,:,index*segLength+1:(index+1)*segLength) = salMapTime(:,:,count0+1:count0+segLength);

%            salMap=salMap(:,:,count0+1:count0+segLength);
%             salMapY=salMapY(:,:,count0+1:count0+segLength);
%             salMapI=salMapI(:,:,count0+1:count0+segLength);
%             salMapQ=salMapQ(:,:,count0+1:count0+segLength);
%             salMapTY=salMapTY(:,:,count0+1:count0+segLength);
%             salMapTI=salMapTI(:,:,count0+1:count0+segLength);
%             salMapTQ=salMapTQ(:,:,count0+1:count0+segLength);
%             salMapSY=salMapSY(:,:,count0+1:count0+segLength);
%             salMapSI=salMapSI(:,:,count0+1:count0+segLength);
%             salMapSQ=salMapSQ(:,:,count0+1:count0+segLength);
        end

        
%         % Display/save saliency maps with original video frames.
%         % Need to change if multiscale in time is incorporated.
%         %figure;
%         for k=1:size(salMap,3)
%             tmpY=imresize(salMapY(:,:,k),[frameH frameW],'bilinear');
%             tmpI=imresize(salMapI(:,:,k),[frameH frameW],'bilinear');
%             tmpQ=imresize(salMapQ(:,:,k),[frameH frameW],'bilinear');
%             tmpA=imresize(salMap(:,:,k),[frameH frameW],'bilinear');
%             tmpYQ=cat(1,tmpY,tmpQ);
%             tmpIA=cat(1,tmpI,tmpA);
%             
%             frameYQ=cat(1,frameY(:,:,count0+k),frameQ(:,:,count0+k));
%             frameIY=cat(1,frameI(:,:,count0+k),frameY(:,:,count0+k));
%             
%             tmpYQIA=cat(2,tmpYQ,tmpIA);
%             frameYQIY=cat(2,frameYQ,frameIY);
%             
%             out=sc(cat(3,tmpYQIA,frameYQIY),'prob_jet');
%             % Display w/o saving
%             %sc(cat(3,tmp,frameY(:,:,count0+k)),'prob_jet');
%             %pause(.03);
%         end
        clear salMap salMapTime;
        clear salMapTY salMapTI salMapTQ;
        clear salMapSY salMapSI salMapSQ;
        index = index+1;
        
        % Reset data matrices and counters, keeping some data for next segment
        if (j~=jStop)
            dataY=dataY(:,:,count-2*hw_t-2*hw_ft+1:count);
            dataI=dataI(:,:,count-2*hw_t-2*hw_ft+1:count);
            dataQ=dataQ(:,:,count-2*hw_t-2*hw_ft+1:count);
            frameY=frameY(:,:,count-2*hw_t-2*hw_ft+1:count);
            count0=hw_t+hw_ft;
            count=2*hw_t+2*hw_ft;
        end
    end
end

%% End of program
