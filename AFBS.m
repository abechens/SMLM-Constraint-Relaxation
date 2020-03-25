function [xopt, infos ] = AFBS(func, paramsResol)
%AfBS is the Nonmonotone Accelerated Proximal gradient algorithm which was
%presented in Accelerated proximal gradient methods for nonconvex
%programming by Li, Huan and Lin, Zhouchen in Advances in neural
%information processing systems 

%   Inputs: 
%   func containts different functions needed:
%       -cost which is the cost function
%       -grad which calculates the gradient 
%       -prox which calculates the proximal
%   paramsResol contains different parameters needed:
%       -xo, the initialization
%       -stopConvX, stopConvX and
%       -itermax are stopping crieterias
%       -eta is the value corresponding to eta in the article
%       -delta is the value corresponding to the delta in the article
%       -norminv if we normalize A. 
%
%   Outputs:
%   xopt: the optimal x
%   infos, which contains the cost development, the number of iterations
%   and the time used.


% -- Set default values

if ~isfield(paramsResol,'itermax'), paramsResol.itermax=8000; end
if ~isfield(paramsResol,'stopConv'),paramsResol.stopConv = 1e-5; end

tic;
infos{1} = [];
infos{2} = [];
infos{3} = [];
% === Définition des variables
etaP = paramsResol.eta;
deltaP = paramsResol.delta;


nonConv=1;
xit = paramsResol.x0;
xold=xit;
zit = xit;
ckit = func.cost(xit);
qkit = 1;
tkold = 1;
tkit = 1;
costFBS= zeros(1, paramsResol.itermax);
it = 0;
proxIndanc=1;

while nonConv
    it = it+1;
    
    yit = xit + tkold/tkit * (zit-xit) + (tkold-1)/tkit * (xit - xold); 
    
    xold = xit;
    % ===== FBS ====
    
    gradientF = func.grad(yit);
    stepsY = yit -  gradientF; %Step is included in gradientF
    
    
    
    % === Prox S^2 ===
    [ zit, proxIndanc ] = func.proximal(stepsY, proxIndanc);
    % === Calcul du cout FBS ====
    
   
    costTemp = func.cost(zit);
    
   % === Else condition ===
   if costTemp <= ckit - deltaP* norm(zit-yit,'fro')^2
       xit = zit;
       costFBS(it) = costTemp;
   else
        gradientF = func.grad(xit);
        stepX = xit -  gradientF;
    
         % === Prox S^2 ===
        [ vkit, proxIndanc ] = func.proximal(stepX, proxIndanc);
        costTempV = func.cost(vkit);
       if costTemp <= costTempV
           xit = zit;
           costFBS(it) = costTemp;
       else
           xit = vkit;
           costFBS(it) = costTempV;
       end
   end
   
    tkold = tkit;
    tkit = (sqrt(4*tkit^2+1 )+1)/2;
    qkit = etaP*qkit +1;
    ckit = (etaP * qkit * ckit+costFBS(it)) / qkit;
    % ==== Convergence FBS ===
   
    if it > paramsResol.itermax
        nonConv = 0;
    elseif it~=1
        diffCout = abs(costFBS(it) - costFBS(it-1))/abs(costFBS(it)+eps);
        diffX = norm(reshape(xit-xold,[],1),'fro')/(norm(xold(:),'fro')+eps);
        if (diffCout < paramsResol.stopConvC && diffX < paramsResol.stopConvX )
            nonConv = 0;
        end
    end
    
end





x=func.failsafe(xit); %Test if converged correctly

if isfield(paramsResol,'norminv') && all(all(x==xit)) %If converged correctly, we reconstruct the proper intensity from the columns of A_i
xopt = paramsResol.norminv.*xit;
else 
    xopt=xit;
end
infos{1} = costFBS(1:it);
infos{2} = it;
infos{3} = toc;
end


