# SMLM-Constraint-Relaxation
Single Molecule localization microscopy code based on a deconvolution and a sparsity parameter L0. Here we use a relaxation of the constrained L0 term, which we note Q(x). 

Main less the user test the file. We minimize the squared norm and the relaxation Q(x), using Nonmonotone Accelerated Proximal gradient algorithm which was presented in Accelerated proximal gradient methods for nonconvex programming by Li, Huan and Lin, Zhouchen in Advances in neural information processing systems. 
In order to do so, we use gradient of the square norm and the proximal of Q(x). There are two functions for the proximal operator, and reccomend the fast one for the calculations, but the normal one to better understand what the program is doing. 
A failsafe function is added to ensure that we always obtain a k-sparse solution.
The file uses a simulated acquistion, Testimage. 
