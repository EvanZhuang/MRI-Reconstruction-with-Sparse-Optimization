# MRI-Reconstruction-with-Sparse-Optimization
Magnetic resonance imaging (MRI) images are known to be sparse. This is an implementation using non-convex penalty function that encourages sparsity.

The penalty function is chosen as the minimax concave penalty (MCP), the algorithm (GIST) can be checked from:

A General Iterative Shrinkage and Thresholding Algorithm for Non-convex Regularized Optimization Problems by Pinghua Gong, Changshui Zhang, Zhaosong Lu, Jianhua Huang, Jieping Ye https://arxiv.org/abs/1303.4434

Run main.m directly and you will see the comparison between popular methods and this implementation.

The Randon transform code and back projection to DFT code are written by Mark Bangert.

![Sample Image](https://github.com/EvanZhuang/MRI-Reconstruction-with-Sparse-Optimization/blob/master/4Compare.jpg)

The solvers are also incuded in the solver folder, select the one you need. GIST_MCP.m used proximal gradient method with Barzilai-Borwein step size, GIST_MCP_Nesterov.m used proximal gradient method with Nesterov acceleration. Remember to put the corresponding subroutine with the solver.

There is detailed explanation of the Nesterov accelerated proximal gradient algorithm with restart that truly guarantees convergence, here:

Linear Convergence of Proximal Gradient Algorithm with Extrapolation for a Class of Nonconvex Nonsmooth Minimization Problems by Bo Wen, Xiaojun Chen, Ting Kei Pong https://arxiv.org/pdf/1512.09302.pdf
