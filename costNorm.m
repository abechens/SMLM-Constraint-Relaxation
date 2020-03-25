function [ cost ] =costNorm(x,paramsConv,Iacq, norminv)
%costNorm calculates the ||Ax-d||^2_2
%   The code is inspired by:
%   https://github.com/esoubies/SMLM-CEL0/blob/master/ComputeNorm_ai.m
%
%   Inputs:
%   x is the x
%   paramsConv contains A, the gaussian kernel, and M, the reduction
%   operator, which functions as this MxM'. 
%   norminv, the matrix contating the inveresed values of the
%   observation matrix
%
%   Returns the value of ||Ax-d||^2_2. 

Afft = fft2(fftshift(paramsConv.A)); % Gaussian Kernel
Aconj = conj(Afft);
Mech = paramsConv.M;
MechT = paramsConv.M';
ech = size(Afft,1)/size(Iacq,1);



xitfft = fft2(norminv.*x);
cost= 0.5*norm(Mech*real(ifft2(Afft.*xitfft))*MechT-Iacq,'fro')^2;

end

