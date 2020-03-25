function [ x ] = failsafe(xin,k,failsafefunc)
%failsafe test if the sparsity is reached. 
%-  If the sparsity is reached, x=xin;
%   if not, an message will appear, xin is projected to the nearest l_0<=k
%   space. After this the optimal values of x is calculated.
%  
%   Inputs:
%   - xin is the vector
%   - k is the sparsisty parameter
%   -failsafefunc is a function that, given a support, calculate the
%   optimal intensity for that support.
%
%   Output:
%   x, which is a minimizer of G_k

if nnz(xin)>k
    disp('===== FAIL-SAFE ACTIVATED======')
    fprintf('===== THE ALGORITHM FAILED TO CONVERGE TO a %d - sparse solution======\n',k)
    fprintf('===== x has sparsity= %d ======\n',nnz(xin))
    xl0 = proxsl0(xin,k);
    x=failsafefunc(xl0);
else
    x=xin;
end

