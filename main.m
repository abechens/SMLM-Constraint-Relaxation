%% Demonstration of algorithm
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%% Parameters to create the observation matrix 
% (We will load it instead of creating it. 

% paramsAcq.ech = 4; % The L factor. How much upsampling
% paramsAcq.N = 64*paramsAcq.ech; %Size of ground truth grid
% paramsAcq.sig = 109.70; % PSF estiamtion (sigma)
% paramsAcq.c = 0; % center of PSF
% paramsAcq.I0 = 1; % intensisty PSF
% paramsAcq.taillePixel = 100/paramsAcq.ech ; % size of each pixel in nm
% 
% %% Creating the Matrix. 
% 
% matriceConv = ParamAcq(paramsAcq);
load('MatriceA.mat') %The matrix A, containing A and M, A is the
%                       convolution matrix, and M ...
%                       the downsampling matrix. (Downsampling done by
%                       M*x*M'). 

%% Defining the parameters to solve to problem
% Converging parameters
paramsResol.stopConvC = 10^(-5);
paramsResol.stopConvX = 10^(-8); %Recommend somewhat higher here than stopConvC as the algorithm may be slow. 
paramsResol.itermax =10^5;

% Parameters connnected to the proximal operator
% IMPORTANT: This is a speed up of the proximal operator, and is not
% important if you choose to use the "slow" proximal operator
paramsResol.sizebetter=150;
paramsResol.alphas = 1;
paramsResol.eta= 0.8; 
paramsResol.delta= 0.1; 

%initialization 
paramsResol.x0 = zeros(size(matriceConv.A) ); %This yields the best results

%Sparsity parameters
paramsResol.kmax = 99; % The number of non-zeros pixels reconstructed

load('testimage.mat') %The test image is the first acquisition of a 361 High
%                       -density acquisition. Simulated, and found at the
%                       page: http://bigwww.epfl.ch/smlm/datasets/index.html

[step, paramsResol.norminv] = calculateParams(matriceConv,paramsResol.alphas);

%% Relaxed algorithm
observation = testimage/max(testimage(:));
func.cost = @(x) costNorm(x,matriceConv,observation,paramsResol.norminv)+costQ(x,paramsResol.kmax)+paramsResol.alphas/2*norm(x(x<0),'fro')^2;
func.grad= @(x) GradientSquare(x,matriceConv,observation,paramsResol.alphas,step,paramsResol.norminv);
func.proximal =@(x,y) proxsQFast(x,step,paramsResol.kmax, y, paramsResol.sizebetter);
func.optimalValue =@(x) optimalValue(x, matriceConv,paramsResol.kmax,observation);
func.failsafe = @(x) failsafe(x,paramsResol.kmax, func.optimalValue);
[GQ, info] = AFBS( func, paramsResol);  
GQ= GQ*max(testimage(:));



