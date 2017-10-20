function [u, v] = optFlowLK(frame1, frame2, window)
% Stolen from:
% https://www.mathworks.com/matlabcentral/fileexchange/ ...
% 48744-lucas-kanade-tutorial-example-1?focused=3854179&tab=example

w = round(window/2);

% Lucas Kanade Here
% for each point, calculate I_x, I_y, I_t
Ix_m = conv2(frame1,[-1 1; -1 1], 'valid'); % partial on x
Iy_m = conv2(frame1, [-1 -1; 1 1], 'valid'); % partial on y
It_m = conv2(frame1, ones(2), 'valid') + conv2(frame2, -ones(2), 'valid'); % partial on t
u = zeros(size(frame1));
v = zeros(size(frame2));

% within window ww * ww
for i = w+1:size(Ix_m,1)-w
   for j = w+1:size(Ix_m,2)-w
      Ix = Ix_m(i-w:i+w, j-w:j+w);
      Iy = Iy_m(i-w:i+w, j-w:j+w);
      It = It_m(i-w:i+w, j-w:j+w);

      Ix = Ix(:);
      Iy = Iy(:);
      b = -It(:); % get b here

      A = [Ix Iy]; % get A here
      nu = pinv(A)*b; % get velocity here

      u(i,j) = nu(1);
      v(i,j) = nu(2);
   end
end

% % downsize u and v
% u_deci = u(1:10:end, 1:10:end);
% v_deci = v(1:10:end, 1:10:end);
% % get coordinate for u and v in the original frame
% [m, n] = size(frame1);
% [X,Y] = meshgrid(1:n, 1:m);
% X_deci = X(1:20:end, 1:20:end);
% Y_deci = Y(1:20:end, 1:20:end);


end

