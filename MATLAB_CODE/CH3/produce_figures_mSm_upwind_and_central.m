%% Script that produces Figures 3.1, 3.2, 3.3, 3.4 and 3.5 from the Thesis:
%
%  Echeverria - Iterative Solution of Convection-Diffusion Problems - 2020
%
%   Written at the Technische Universitaet Berlin.
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
%   a subfolder named 'figures' in the working directory.
%
%  subordinate functions:
%          
%        shishkin_experiments_upwind.m
%        shishkin_experiments_central.m
%
% Written by Carlos Echeverria on October 10, 2019.
% Edited  by C.E. on February 20, 2020.


   alpha = 1.0;  % scalar convection coefficient
    beta = 0.0;  % scalar reaction coefficient
       N = 66;  % number of intervals to subdivide the domain [0,1]
  scaled = 0;    % choose 1 for using scaled versions of coeff. matrix

% epsilons = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]; % desired perturbation parameters
epsilons = [1e-8,1e-6,1e-4]; % desired perturbation parameters


for j = 1:length(epsilons)

    epsilon = epsilons(j);

    disc = 1;
    clf;
    shishkin_experiments_upwind(epsilon, alpha, beta, N, disc, scaled);
    mkdir figures
    file_name = sprintf('/figures/mSm_upwind_eps_%5.0e_N_%d.eps', epsilon, N);
    saveas(gcf,[pwd file_name],'epsc')
    
    disc = 2;
    clf;
    shishkin_experiments_central(epsilon, alpha, beta, N, disc, scaled);
    file_name = sprintf('/figures/mSm_central_eps_%5.0e_N_%d.eps', epsilon, N);
    saveas(gcf,[pwd file_name],'epsc')
    

end
