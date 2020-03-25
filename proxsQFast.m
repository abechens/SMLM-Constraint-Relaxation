function [ prox4, rightInd ] = proxsQFast(xin,step,k, proxIndanc, sizebetter)
%Calculate the proximal operator of the relaxation Q(x) 
%   This is the fast  version.
%- xin is the vector
%- step is the step size. 
%- k is the sparsisty parameter
%- proxIndanc is a parameter which decides where to start looking for the
%\tau. We can suppose it is around the same intervall as the previous
%iteration.
%- sizebetter is the number of indexes we test each time. 

rhot =1/step;
xin2 = xin(:);

[xsort,xind] = sort(abs(xin2),'descend');
prox3=zeros(size(xin2));
xsort2=xin(xind); %xsort2 is x vectorized and sorted by its magnitude
len=length(xin2);

%% for loop to search optimal \tau if \rho* y_k+1>y_k
tausug=xsort(k);
rightInd=k;
if xsort(k) <= xsort(k+1) * rhot % Check if we need to search
    xprem = xsort(xsort(1:k) <=  xsort(k+1) * rhot); % All entries smaller than the largest on the 'right side'
    xfin= rhot * xsort([false(k,1); rhot * xsort(k+1:end)>=xsort(k)]); % All entries larger than the smallest left side (draw)
  %%Keeping only uniques values
    [xpremN, xidxu, xidxc] = unique(xprem,'stable'); %keeping only unique values.
    [xpremcount, ~, ~] = histcounts(xidxc,numel(xidxu)); %saving the number of occurency of each value
    [xfinN, xidxu, xidxc] = unique(xfin,'stable');%keeping only unique values.
    [xfincount, ~, ~] = histcounts(xidxc,numel(xidxu));%saving the number of occurency of each value
    yunique = sort([xpremN; xfinN],'descend'); %sort all values 
    indices= 1: length( yunique ) -1;
   
    testLength= length(indices);
    proxSplit.large=min(proxIndanc,testLength); %not to overstep
    proxSplit.small=min(proxIndanc,testLength);%not to overstep
    if testLength==1 %only one possiblity
        min1 = yunique(indices+1); %
        max1 = yunique( indices);% Searching tau in the intervall [min1, max1]
        [ n1, n2, xvecT, xvecB ] = countUnderOver(min1, max1, xpremN, xpremcount,xfinN,xfincount );
        
        tausug = (rhot * ( xvecT) +  (xvecB))...
            / ( rhot * n1 + n2); %XvecB already multiplied by rhot. 
        rightInd =indices;
    else %this means we have to search intervals
        notConv=1;
        direction=-1; %We assume that \tau is in a group of indeces less then last
        while notConv
            [ indicenew, proxSplit ] = indicetest(indices,proxSplit, sizebetter, direction);
            %Deciding which intervall to search in, and saving that where
            %we search
            
            for ftau = indicenew
                min1 = yunique(ftau+1);
                max1 = yunique( ftau);
                [ n1, n2, xvecT, xvecB ] = countUnderOver(min1, max1, xpremN, xpremcount,xfinN,xfincount );
                
                tausug = (rhot * ( xvecT) +  (xvecB))...
                    / ( rhot * n1 + n2);
                
                if tausug+eps>=min1 && tausug-eps<=max1 % Condition satisfied
                    notConv=0;
                    rightInd =ftau;
                    break
                end
            end
            if proxSplit.small~=1 && proxSplit.large~=length(yunique)-1
                direction = -1*direction; %looking in the oposite direction of the start, we  alternate between seaching over our under
            else
                if proxSplit.small==1 %If we cannot go any lower
                    direction=1;
                else
                    direction=-1;
                end
            end
            
        end
    end
end


%Finding the index, that is k* and k**
index1 = [abs(xsort2(1:k)) >= tausug;false(len-k,1)]; % Stay unchanged
index2 = [false(k,1);rhot * abs(xsort2(k+1:end))<=tausug]; % those whos value is 0
index3 = logical(true(len,1)-index1 - index2); % the others

prox3(index1) = xsort2(index1);
prox3(index2) = 0;
prox3(index3) = (rhot*xsort2(index3)-sign(xsort2(index3)).*tausug)/(rhot - 1);


prox4(xind) = prox3; %Reordering back

prox4 = reshape(prox4,size(xin)); %Reshaping


end

