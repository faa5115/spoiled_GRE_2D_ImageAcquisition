function res = ifftnc(x)
%
%
% res = ifftnc(x)
% 
% orthonormal centered N-dimensional ifft
%

res = ifftshift(ifftn(ifftshift(x)));
res = sqrt(length(x(:)))*res;

