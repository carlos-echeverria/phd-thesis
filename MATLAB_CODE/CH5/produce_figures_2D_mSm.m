% Script for generating the figures showing the convergence behavior of 
% the multiplicative Schwarz method as well as the theoretical convergence
% bounds presented in Chapter 5  (Figures 5.1 and 5.2)of the PH.D. thesis:
%
%  [ Echeverria - Iterative solution of discretized convection-diffusion 
%                 problems - Technische Universitaet Berlin - 2020 ]
%
%   This script perfroms numerical experiments for solving the system:
%   Ax = b, where the matrix A is obtained from the discretization of the
%   following 2D conv-diff problem posed on a Shishkin mesh with one
%   baoundary layer:
%
%  -eps*Delta(u) + (0,1).Nabla(u) + u = f in Omega,  u = g on Gamma,
%
%   The system is solved via the multilpicative Schwarz method using the 
%   function: experiments_mSm_upwind.m.
%   For each set of parameters epsilon, N, and M, the function produces a 
%   Figure that plots the error produced by the method at each iteration 
%   step. The theoretical bounds presented in the thesis are also plotted.
%   The figures are then stored in an eps file in a subfolder of the
%   current working directory.
%
% Written by Carlos Echeverria on October 10, 2019.
% Edited  by C.E. on February 25, 2020.


prob = 1;       %  1: Conv-Diff, 2: Poisson
   N = 30;      %  number of intervals of the grid in x-dir (must be even)
   M = 40;      %  number of intervals of the grid in y-dir (must be even)

epsilons = [1e-8, 1e-6, 1e-4, 1e-2]; % desired perturbation parameters

for j = 1:length(epsilons)

    epsi = epsilons(j);

    [rho, rho_e, max_err] = experiments_mSm_upwind(prob, N, M, epsi);
    max_err_vec(j) = max_err;
    mkdir figures
    file_name = sprintf('/figures/mSm_upwind2D_eps_%5.0e_N_%d_M_%d.eps', epsi, N, M);
    saveas(gcf, [pwd file_name],'epsc')

end

%% output values to console

format shorte; matrix_size = (N-1)*(M-1); epsilons, max_err_vec

