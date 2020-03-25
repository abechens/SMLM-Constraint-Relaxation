function [ S ] = costQ( xin,k )
%costQ Calculate the value Q(x).
%
%   Inputs:
%   xin is x
%   k is the sparsity parameter.
%
%   Output:
%   S, the value of Q(x)

%Vectorize and sort by the absolute value

xin2 = xin(:);
[xsort] = sort(abs(xin2),'descend');

%%%% Search for T(x) %%%% 

if k==1 %Easy case since 1<=T(x)<=k
    tcor=1;
else
    
    for tsug=1:k
        if tsug==k %If all other option are tested
            
            break;
            
        else
            %Test of condition, added with eps since x tends to be small
            %which add small inprecisions.
            if 1/tsug * sum(xsort(k-tsug+1:end))<=xsort(k-tsug)+eps && 1/tsug * sum(xsort(k-tsug+1:end))+eps>=xsort(k-tsug+1)
                
                break;
            end
        end
    end
    tcor = tsug; %T(x) is found. 
end
S=-1/2 * sum(xsort(k-tcor+1:end).^2)+1/(2*tcor) *   sum(xsort(k-tcor+1:end))^2; 
end