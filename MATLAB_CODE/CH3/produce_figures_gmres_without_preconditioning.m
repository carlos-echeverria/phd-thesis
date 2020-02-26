%% Script that produces Figures 2.8 and 2.9 from the Thesis:
%
%  Echeverria - Iterative Solution of Convection-Diffusion Problems - 2020
%
%   Written at the Technische Universitaet Berlin.
%
%   This script perfroms numerical experiments for solving the system:
%   Ax=b, where the matrix A is obtained from the discretization of the
%   following 1D conv-diff problem posed on a Shishkin mesh with one
%   transition point:
%
%   -epsilon * u'' + alpha * u' + beta * u = f in [0,1],  u(0) = u(1) = 0.
%
%   The system is solved via GMRES without preconditioning using the
%   function: gmres_shishkin.m. For each set of parameters: epsilon, alpha,
%   beta, and N, the function produces a Figure that compares the
%   performance of GMRES when the matrix is obtained via the upwind
%   finite diference scheme as well as the central diference scheme.
%   Performance for the scaled versions of the resulting matrices is also
%   included in the comparison. The figure is then saved in an eps file in
%   the working directory.
%
%  subordinate functions:
%          
%        gmres_shishkin.m
%
% Written by Carlos Echeverria on October 10, 2017.
% Edited by C.E. on February 20, 2020.



alpha = 1.0;  % scalar convection coefficient
 beta = 0.0;  % scalar reaction coefficient
    N = 198;  % number of intervals to subdivide the domain [0,1]

epsilons   = [1e-8, 1e-6, 1e-4, 1e-2]; % desired perturbation parameters

for j = 1:length(epsilons)

    epsilon = epsilons(j);

    gmres_shishkin(epsilon,alpha,beta,N);
    mkdir figures
    file_name = sprintf('/figures/gmres_eps_%5.0e_N_%d.eps',epsilon, N);
    saveas(gcf,[pwd file_name],'epsc')
    clf;

end
