% function v = convolve(u, k)
% computes 2D convolution u * k via fft
% k is PSF / OTF in Fourier space
% last modified : 26.11.2015
%

function v = convolve(u, k)
if nargin < 2
   error('Not enough arguments!');
end
sk=size(k);
assert(numel(u)==numel(k),'dimension mismatch!');

if all(size(squeeze(u))==size(k))
   v=ifftshift(real(ifft2(fft2(squeeze(u),sk(1),sk(2)) .* k)));
else
   w=reshape(u,size(k));
   v=ifftshift(real(ifft2(fft2(w) .* k)));
   v=reshape(v,size(u));
end

end

