function [ Grad ] =GradientSquare(x,paramsConv,Iacq,alphas,step, norminv)
%GradientSquare calculates the gradient of the square norm + the gradient
%of the distance function from R+. 
%
%   Inputs:
%   x
%   paramsConv, which contains the gaussian kernal, A, and the reduction of
%   dimension
%   Iacq the observed data
%   step, the gradient step length
%   norminv, the matrix contating the inveresed values of the
%   observation matrix
%   
%   Outputs:
%   Grad, the gradient



Afft = fft2(fftshift(paramsConv.A)); 
Aconj = conj(Afft);

xfft = fft2(norminv.*x);
xtemp=x;
xtemp(x>0)=0;
Grad =step.*( norminv.*real(ifft2(Aconj.*fft2((paramsConv.M'*(paramsConv.M * real(ifft2(Afft .* xfft)) *paramsConv.M' - Iacq)*paramsConv.M))))+alphas*xtemp);

end

