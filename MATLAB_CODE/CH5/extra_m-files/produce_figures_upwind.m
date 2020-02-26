%% Script that produces Figures 3.1, 3.2, 3.3, 3.4 and 3.5 from the Thesis.
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
%   the working directory.
%
% Written by Carlos Echeverria on October 10, 2019.

prob = 1;       %  1: Conv-Diff, 2: Poisson
   N = 30;      %  number of intervals of the grid in x-dir (must be even)
   M = 40;      %  number of intervals of the grid in y-dir (must be even)


epsilons = [1e-8, 1e-6, 1e-4, 1e-2]; % desired perturbation parameters

for j = 1:length(epsilons)

    epsi = epsilons(j);

    % % uncomment for MSM experiments
    [rho, rho_e, max_err] = experiments_mSm_upwind(prob, N, M, epsi);
    max_err_vec(j)=max_err;
    file_name = sprintf('mSm_upwind2D_eps_%5.0e_N_%d_M_%d.eps', epsi, N, M);
    saveas(gcf, file_name,'epsc')

    % % uncomment for preconditioned GMRES experiments
%     experiments_gmres_precond_upwind(prob, N, M, epsi);
%     file_name = sprintf('gmres_precond_upwind2D_eps_%5.0e_N_%d_M_%d.eps',epsi, N, M);
%     saveas(gcf,file_name,'epsc')

    % % uncomment for unpreconditioned GMRES experiments
%     experiments_gmres_upwind(prob, N, M, epsi);
%     file_name = sprintf('gmres_upwind2D_eps_%5.0e_N_%d_M_%d.eps',epsi, N, M);
%     saveas(gcf,file_name,'epsc')


end

format shorte;
problem_size=(N-1)*(M-1)
epsilons
max_err_vec
