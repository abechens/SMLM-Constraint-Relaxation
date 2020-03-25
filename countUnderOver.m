function [ n1, n2, xvecT, xvecB ] = countUnderOver(min1, max1, xpremN, xpremcount,xfinN,xfincount )
%countUnderOver yields the sum of all y_i in n_1 and y_i in n_2, as well as n_1 and n_2.
%   The code is used to find \tau in the proximal operator
%
%   Inputs:
%   min1 is the lower intervall in which we search \tau
%   max1  is the upper intervall in which we search \tau
%   xpremN is the sorted unique vector of the of the k first
%   xpremcount is the number of each element 
%   xfinN is the sorted unique vector of the of the for i>k
%   xfincount is the number of each eleemtn
%
%   Output:
%   n1: the size of n_1
%   n2: the size of n_2
%   xvecT: is the sum of |x_i^\downarrow| in the interval n_1
%   xvecB: is the sum of |x_i^\downarrow| in the interval n_2


    indn1= xpremN<=min1;
    n1 = sum(xpremcount(indn1));
    indn2 = xfinN >= max1 ;
    n2 = sum(xfincount(indn2));
    
    xvecT = sum(xpremN(indn1).*xpremcount(indn1)');
    xvecB = sum(xfinN(indn2).*xfincount(indn2)');

end

