function res = fftnc(x)

% res = fftnc(x)
% 
% orthonormal forward N-dimensional FFT
%

res = fftshift(fftn(fftshift(x)));
res = res/sqrt(length(x(:))); % normalization

