% Script for generating Table 5.1 presented in Chapter 5 
% of the PH.D. thesis:
%
%  [ Echeverria - Iterative solution of discretized convection-diffusion 
%                 problems - Technische Universitaet Berlin - 2020 ]
%
%   This script obtains the coefficient matrices from the discretization 
%   of the following 2D conv-diff problem posed on a Shishkin mesh with 
%   one boundary layer:
%
%  -eps*Delta(u) + (0,1).Nabla(u) + u = f in Omega,  u = g on Gamma,
%
%   The matrices are then used to calculate the "exact" as well as the
%   theoretical bounds for the convergence factor of the multiplicative
%   Schwarz method using the function: experiments_mSm_upwind.m 
%   For each set of parameters epsilon, N, and M, the function produces a 
%   Figure that plot the relative error produced by the method at each
%   iteration step (w.r.t. the backsalsh solution). The figure is then 
%   stored in an eps file in a subfolder of the current working directory
%   and the values of the table are shown in the console.
%
%
% Written by Carlos Echeverria on October 10, 2019.
% Edited  by C.E. on February 25, 2020.

clearvars;
    prob = 1;     %  1: Conv-Diff, 2: Poisson
epsilons = [1e-8, 1e-6, 1e-4, 1e-2];
%  sizes = 2.^[3:1:6]+2;
  sizesx = [20,20,30,30,40,40,50,50];
  sizesy = [20,30,30,40,40,50,50,60];
% sizesx = [4];
% sizesy = [8];

for i=1:length(sizesx)
    
    N = sizesx(i);
    M = sizesy(i);
    
    for j = 1:length(epsilons)

        epsi = epsilons(j);
        [rho_bnd, rho_e_inf, max_err] = experiments_mSm_upwind(prob, N, M, epsi);
        
        mkdir figures
        file_name = sprintf('/figures/mSm_upwind2D_eps_%5.0e_N_%d_M_%d.eps', epsi, N, M);
        saveas(gcf, [pwd file_name],'epsc');

        fprintf('N:%d, M:%d, A-size:%dx%d, epsi:%5.1e, rho_e:%5.1e, bound:%5.1e. \n', N,M,(N-1)*(M-1),(N-1)*(M-1),epsi, rho_e_inf, rho_bnd);

    end

end
