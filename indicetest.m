function [ indicenew, proxSplit ] = indicetest(indices,proxSplit, sizebetter, direction)
%indicetest gives the a groupe of new indices to test and updates
%proxIndanc which contains the the lower and upper bounds of indices tested
%
%   Inputs:
%   Indices is the set of indices to check, 
%   proxIndanc containts borns tested
%   sizebetter is each "jump"
%   direction is if we are looking it to go down or up in the index.
%
%   Output: 
%   indicenew, the groupe of indices to test if \tau is in it
%   proxIndanc, containing 
%               -large, which is the upper bound of indices already tested
%               -small, wich is the lower bound of indices already tested
if direction ==1
    maxdist=min(length(indices),proxSplit.large+sizebetter);
    indicenew= indices(proxSplit.large:direction:maxdist);
    proxSplit.large=maxdist+1;
else
    maxdist=max(1,proxSplit.small-sizebetter);
    indicenew= indices(proxSplit.small:direction:maxdist);
    proxSplit.small=maxdist;
end

end

