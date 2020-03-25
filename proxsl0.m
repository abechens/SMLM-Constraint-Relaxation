function [ prox3 ] = proxsl0(xin,k)
% Proximal of the hard threshold
%   Keeps only the k biggest values
%   Inputs:
%   x, the vector on which we perform the proximal operator on
%   k, the sparsity constraint.
%
%   Output:
%   prox3, the proximal operator

xin2 = xin (:);
[~,xind] = sort(abs(xin2),'descend'); %Sort the vector by its absolute value
xsort=xin2(xind); %sort them  but keeping the sign
prox2=zeros(size(xin2)); %This will be the proximal opertor sorted

prox2( 1:k ) =xsort( 1:k ); %Keeping the k largest value
prox2(k+1:end) = 0;



prox3(xind) = prox2; % sort back to right position

prox3 = reshape(prox3,size(xin)); %and reshape.


end

