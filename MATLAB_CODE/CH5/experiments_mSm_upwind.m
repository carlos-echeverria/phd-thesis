function [rho_bnd, rho_e_inf, max_err] = experiments_mSm_upwind(prob, N, M, epsi)
%   This script perfroms numerical experiments for solving the system:
%   Ax = b, where the matrix A is obtained from the discretization of the
%   following 2D conv-diff problem posed on a Shishkin mesh with one
%   baoundary layer near the outflow boundary:
%
%  -eps*Delta(u) + (0,1).Nabla(u) + u = f in Omega,  u = g on Gamma,
%
%   The system is solved via the multilpicative Schwarz method.
%   For each set of parameters epsilon, N, and M, the function produces a 
%   Figure that plots the error produced by the method at each iteration 
%   step. The theoretical bounds presented in the thesis are also plotted.
%
%   function call: 
%
% [rho_bnd, rho_e_inf, max_err] = experiments_mSm_upwind(prob, N, M, epsi)
%
%   input:
%        prob: type of problem. 1 for conv-diff, 2 for Poisson.
%           N: number of intervals of the Shishkin mesh in the x-direction.
%           M: number of intervals of the Shishkin mesh in the y-direction.
%        epsi: scalar diffusion coefficient (perturbation parameter).
%
%   output:
%       rho_bnd: theoretical bound from the Thesis.
%     rho_e_inf: theoretical convergence factor using inf-norm
%       max_err: maximum discrepancy between exact (backslash) and 
%                computed solution (mSm).
%
%
% Written by Carlos Echeverria on October 10, 2019.
% Edited  by C.E. on February 25, 2020.

%% Define specific problem (parameters must be hardcoded in function below)

[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);

%% Create Shsihkin mesh
  
    [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );

%% Get coefficien matrix of linear system
  
  [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);

%% Get right hand side vector of linear system

  [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
  
%% Solve system using MATLAB's backslash operator

[SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );

%% Solve system using the multiplicative Schwarz method

SOLVER.order = 1; % order 1: T=(I-P2)(I-P1)
[SOLUTION] = block_Multiplicative_Schwarz(SOLVER, MESH, SOLUTION);

% creates new fields in SOLUTION structure for storing computed solution
% and residual of ordering T=Q2Q1
if prob == 1
    SOLUTION.err_s1 = SOLUTION.err_s_inf;
elseif prob == 2
    SOLUTION.err_s1 = SOLUTION.err_s_2;
end
  SOLUTION.x_s1   = SOLUTION.x_s;
  SOLUTION.res_s1 = SOLUTION.res_s;

% order 2: T=(I-P1)(I-P2)
SOLVER.order = 2;
[SOLUTION] = block_Multiplicative_Schwarz(SOLVER, MESH, SOLUTION);

% creates new fields in SOLUTION structure for storing computed solution
% and residual of ordering T=Q1Q2
if prob == 1
    SOLUTION.err_s2 = SOLUTION.err_s_inf; % error in inf-norm for conv-diff
elseif prob == 2
    SOLUTION.err_s2 = SOLUTION.err_s_2;   % error in 2-norm  for Laplace
end
  SOLUTION.x_s2   = SOLUTION.x_s;
  SOLUTION.res_s2 = SOLUTION.res_s;


%% Compute Error bounds
[SOLUTION] = block_Schwarz_error_bounds_model_problems(EQ, MESH, SOLUTION);

         rho_bnd = SOLUTION.rho_bnd;
%% Compute exact rho from Lemma 3.3

% [AH,Ah,Z,ZH,Zh,blkA,blkAh,blkAH,blkZ,blkZH,blkZh,Zdev,ZHdev,Zhdev]  = block_Schwarz_extract_and_partition(SOLVER.A, MESH)
[~,~,~,~,~,~,~,~,~,rho_e_inf] = block_Schwarz_compute_rho_e( SOLVER.A, MESH);

%% Plot Results

%  Algebraic error generated by the Schwarz method
block_Schwarz_ploterr(EQ, SOLUTION);

%  Solutions obtained by the different approaches
[max_err] = block_Schwarz_plotsol(MESH, SOLUTION);
figure(202),close;
figure(203),close;
figure(204),close;

end
