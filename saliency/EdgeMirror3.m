function z = EdgeMirror3(x, width)
%
% Code from the Matlab toolbox avialable at 
% https://users.soe.ucsc.edu/~milanfar/research/rokaf/.html/SaliencyDetection.html#Matlab
% which implements the algorithms in the following paper:
% Hae Jong Seo, and Peyman Milanfar, "Static and Space-time Visual Saliency Detection by Self-Resemblance", 
% The Journal of Vision 9(12):15, 1-27, http://journalofvision.org/9/12/15/, doi:10.1167/9.12.15
%
y = cat(2, x(:, width(2)+1:-1:2,:), x, x(: ,end-1:-1:end-width(2),:));
y = cat(1, y(width(1)+1:-1:2, :,:), y, y(end-1:-1:end-width(1), :,:));
z = cat(3, y(:,:,width(3)+1:-1:2), y, y(:,:,end-1:-1:end-width(3)));

