%% Script that produces Table 3.2 from the Thesis.
%
%   This script perfroms numerical experiments for solving the system:
%   Ax = b, where the matrix A is obtained from the discretization of the
%   following 1D conv-diff problem posed on a Shishkin mesh with one
%   transition point:
%
%   -epsilon * u'' + alpha * u' + beta * u = f in (0,1),  u(0) = u(1) = 0.
%
%   The system is solved via GMRES preconditioned with the multilpicative
%   Schwarz method using the function: gmres_shishkin_precond.m.
%   For each set of parameters epsilon, alpha,
%   beta, and N, the function produces a Figure that compares the
%   performance of GMRES when the matrix is obtained via the upwind
%   finite diference scheme as well as the central diference scheme.
%   Performance for the scaled versions of the resulting matrices is also
%   included in the comparison. The figure is then saved in an eps file in
%   the working directory. An avi video file is also saved.
%
% Written by Carlos Echeverria on October 10, 2019.2016
% modified for the thesis by Carlos Echeverria on October 16, 2019

prob = 1;       %  1: Conv-Diff, 2: Poisson
   N = 30;      %  number of intervals of the grid in x-dir (must be even)
   M = 40;      %  number of intervals of the grid in y-dir (must be even)

epsilons   = [1e-8, 1e-6, 1e-4, 1e-2];
local_tols = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
% local_tols = [1e-1,1e-2,1e-3,1e-4, 1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14];

 for i=1:length(epsilons)

     epsilon= epsilons(i);

     for j = 1:length(local_tols)

         local_tol = local_tols(j);

         [res1,res2]= experiments_precon_inexact_table(prob,N,M,epsi,local_tol);
         fprintf('epsi:%5.1e, tol:%5.1e, upwind exact:%d, upwind inexact:%d, \n', epsilon, local_tol, length(res1)-1, length(res2)-1);

     end

 end
