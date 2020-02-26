% Script for generating Tables 5.2, 5.3 and 5.4 as well as Figure 5.7 
% presented in Chapter 5 of the PH.D. thesis:
%
%  [ Echeverria - Iterative solution of discretized convection-diffusion 
%                 problems - Technische Universitaet Berlin - 2020 ]
%
%   This script perfroms numerical experiments for solving the system:
%   Ax = b, where the matrix A is obtained from the discretization of the
%   following 1D conv-diff problem posed on a Shishkin mesh with one
%   transition point:
%
%   -epsilon * u'' + alpha * u' + beta * u = f in (0,1),  u(0) = u(1) = 0.
%
%   The system is solved via GMRES preconditioned with the *inexact* 
%   multilpicative Schwarz method using the function: 
%   experiments_gmres_precond_inexact_upwind.m.
%   For each set of parameters epsi, N, M and local_tol, the function 
%   produces a Figure that compares the performance of the preconditioned 
%   GMRES with exact and inexact local solves. A figure comparing the 
%   obtained solutions using the backslash operator and the preconditioned
%   method with exact and inexact local solves is also produced. The 
%   figures are then stored in eps files in a subfolder of the current 
%   working directory.The values of the tables are printed in the console.
%
%
% Written by Carlos Echeverria on October 10, 2019.
% Edited  by C.E. on February 25, 2020.


prob = 1;       %  1: Conv-Diff, 2: Poisson
   N = 30;      %  number of intervals of the grid in x-dir (must be even)
   M = 40;      %  number of intervals of the grid in y-dir (must be even)

    epsilons = [1e-8, 1e-6, 1e-4, 1e-2];
  local_tols = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];

for i=1:length(epsilons)
         epsi = epsilons(i);
    % vid = VideoWriter('inexact_M30_N40_e-8.avi');
    % open(vid);
    for j = 1:length(local_tols)
        local_tol = local_tols(j);

        experiments_gmres_precond_inexact_upwind(prob, N, M, epsi, local_tol);
        
        mkdir figures
        file_name = sprintf('/figures/gmres_precond_upwind2D_eps_%5.0e_inexact-%5.0e_N_%d_M_%d.eps', epsi, local_tol, N, M);
        saveas(gcf,[pwd file_name],'epsc')
        % frame = getframe(gcf);
        % writeVideo(vid,frame);
    end
% close(vid);
end
