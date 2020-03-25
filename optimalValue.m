function [ xopt ] = optimalValue(x, paramsConv,kmax,Iacq)
%optimalValue creates A_s, the matrix A with only the columns defined by
%the support of x. We calculate xopt by performin x=(A_s'*A_s)\A_s'*Iacq.
%   Input:
%   x, 
%   paramsConv, which contains the gaussian kernel and reduction operator
%   kmax, the sparsity constraint.
%   Iacq, the observation
%
%   Output:
%   xopt


%% Creating A_s
[smalldim,largedim]=size(paramsConv.M);

A=zeros(smalldim^2,kmax);

index = x~=0;
indT = find(index);
for i = 1:kmax
    Ibase = zeros(largedim);
    Ibase(indT(i))=1;
    
    
    
    
    %% Convolution
    
    % ==== Convolution =====
    
    Ibasefft = fft2(Ibase);
    
    Iconv = real(ifft2(Ibasefft.*fft2(fftshift(paramsConv.A))));
    
    % ===== Sreduction =====
    
    IconvSousEch = paramsConv.M*Iconv*paramsConv.M';
    A(:,i)=IconvSousEch(:);
end
% Optimal x_values for the support
xtemp= (A'*A)\A'*Iacq(:);
%reshaping
xopt=zeros(size(x));
xopt(index)=xtemp;
xopt = reshape(xopt,size(x));
end

