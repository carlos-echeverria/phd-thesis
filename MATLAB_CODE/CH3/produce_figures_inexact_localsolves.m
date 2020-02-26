
%
%   ... produce figures of inexact localsolves in preconditioning
%
%% Script that produces Figures 3.8 and 3.59 from the Thesis.
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
% Written by Carlos Echeverria on October 10, 2019.


   alpha = 1.0;  % scalar convection coefficient
    beta = 0.0;  % scalar reaction coefficient
       N = 10;  % number of intervals to subdivide the domain [0,1]
  scaled = 0;    % choose 1 for using scaled versions of coeff. matrix
   order = 1;
    disc = 1;

  epsilons = [1e-6]; % desired perturbation parameters
% local_tols = [1e-1,1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10];
  local_tols = [1e-1,1e-2];

% vid = VideoWriter('inexact_N198_e8.avi');
% open(vid);
for i = 1:length(epsilons)

    epsilon = epsilons(i); % desired perturbation parameters
    
    gmres_shishkin_precon_inexact(epsilon,alpha,beta,N, order, local_tols);
    
    mkdir figures
    file_name = sprintf('/figures/gmres_8_198_inexac.eps');
    saveas(gcf,[pwd file_name],'epsc')
%     frame = getframe(gcf);
%     writeVideo(vid,frame);
end
% close(vid);
