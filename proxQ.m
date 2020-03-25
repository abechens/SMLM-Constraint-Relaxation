function [ prox4 ] = proxQ(xin,step,k)
%Calculate the proximal operator of the relaxation Q(x) 
%  
%- xin is the vector
%- step is the step size. 
%- k is the sparsisty parameter
rhot=1/step;
xin2 = xin(:);

[xsort,xind] = sort(abs(xin2),'descend');
prox3=zeros(size(xin2));
xsort2=xin(xind);
len=length(xin2);

%% for loop to search optimal \tau if \rho/gamm y_k+1>y_k
tausug=xsort(k);

if xsort(k) <= xsort(k+1) * rhot % Check if we need to search
    xprem = xsort(xsort(1:k) <=  xsort(k+1) * rhot); % All entries smaller than the largest on the 'right side'
    xfin= rhot * xsort([false(k,1); rhot * xsort(k+1:end)>=xsort(k)]); % All entries larger than the smallest left side (draw)
%     ynewsort = sort([xprem; xfin],'descend'); % We only need to search in this intervall

    [xpremN, xidxu, xidxc] = unique(xprem,'stable');
    [xpremcount, ~, ~] = histcounts(xidxc,numel(xidxu));
   
   
    [xfinN, xidxu, xidxc] = unique(xfin,'stable');
    [xfincount, ~, ~] = histcounts(xidxc,numel(xidxu));
    yunique = sort([xpremN; xfinN],'descend');

    indic= 1: length( yunique ) -1;
    for ftau = indic
        min1 = yunique(ftau+1);
        max1 = yunique( ftau);
        [ n1, n2, xvecT, xvecB ] = countUnderOver(min1, max1, xpremN, xpremcount,xfinN,xfincount );
       
        tausug = (rhot * ( xvecT) + gamm * (xvecB))...
            / ( rhot * n1 +gamm * n2);
        %         costi(ftau)=costTau(xin,gamm,rhot,k, tausug(ftau));
        %         end
        if tausug+eps>=min1 && tausug-eps<=max1 % Condition satisfied
            break
        end
    end
end


%Finding the index, that is k* and k**
index1 = [abs(xsort2(1:k)) >= tausug;false(len-k,1)]; % Stay unchanged
index2 = [false(k,1);rhot/gamm * abs(xsort2(k+1:end))<=tausug]; % those whos value is 0
index3 = logical(true(len,1)-index1 - index2); % the others

prox3(index1) = xsort2(index1);
prox3(index2) = 0;
prox3(index3) = (rhot*xsort2(index3)-gamm*sign(xsort2(index3)).*tausug)/(rhot - gamm);


prox4(xind) = prox3; %Reordering back

prox4 = reshape(prox4,size(xin)); %Reshaping


end


