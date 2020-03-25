function [ step, norminv ] =calculateParams(paramsConv,alphas)
%calculateParams calculates two important parameters for the code
%
%   The code is inspired by:
%   https://github.com/esoubies/SMLM-CEL0/blob/master/ComputeNorm_ai.m
%
%   Input:
%   paramsConv contains A, the gaussian kernel, and M, the reduction
%   operator, which functions as this MxM'. 
%   alphas, the regularization parameter linked to the positivity
%   
%   Output:
%   step: for the gradient steps, 1 over the lipschitzconstant.
%   norminv, the matrix contating the inveresed values of the
%   observation matrix

Afft = fft2(fftshift(paramsConv.A)); % Noyau gaussien dans Fourier
Aconj = conj(Afft);
Mech = paramsConv.M;
MechT = paramsConv.M';
ech = size(Afft,1)/size(Mech,1);

% === Calcul norme une colonne ===
for i = 1:ech
    for j = 1:ech
        matr = Mech * circshift(paramsConv.A, [j-1,i-1]) *MechT;
        normai(i,j) = sqrt(sum(matr(:).^2));
    end
end

normai = repmat(normai, size(Mech,1));
norminv = 1./normai;

normA= max(norm(Mech*MechT)*abs(Afft(:)))*max(norminv(:));
lipCst =normA^2 +alphas+ 10^-5;
step = 1/lipCst;

end

