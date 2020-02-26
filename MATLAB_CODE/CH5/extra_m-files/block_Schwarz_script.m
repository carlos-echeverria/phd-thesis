% Script for computing the error norms as well as the theoretical bouds
% for the multiplicative Schwarz method applied to a finite difference
% discretization of the singularly perturbed elliptic BVP:
%
% -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma,
%
% where Gamma is the boundary of Omega = [0,1]x[0,1].
%
% In the context of singularly perturbed PDE's, 'eps' is known as the
% scalar diffusion coefficient, 'alpha' = (wx, wy) is the 'wind vector',
% and 'beta' is the scalar reaction coefficient.
%
% The discretization is performed using upwind finite difference operators
% posed on a two-dimensional Shishkin, regular, or mixed mesh.
%
% This script is based on 4 data structures that contain all the values
% needed to describe, discretize, solve, and alayze the BVP; mainly:
% EQ, MESH, SOLVER and SOLUTION.
%
% The values stored in these 4 structures change as the script progresses.
% The user must hardcode the initial values of the parameters for a
% desired specific problem in the m-file: block_Schwarz_getparams.m
% The documentation of that file has a description of the values found
% in each of the data structures.
%
% There are two sets of predefined parameters: one set for the
% convection-diffusion probelm and another set for the Laplace problem.
% Each set of parametes can be chosen by setting the variable 'prob' to
% 1 for the conv-diff problem and 2 for the Laplace problem.
%
% For further reference regarding the theory behind the bounds for the
% error or regarding the discretization procedure see the manuscript:
%
% [ Echeverria, Liesen, Tichy - Analysis of the multiplicative Schwarz
%   method for matrices with a special block structure - 2019 ]
%   manuscript version: 06.09.19
%
% Written by Carlos Echeverria on January 15, 2019.
% Last Edited by C. E. on September 24, 2019.

clearvars; format longe;


%% Obtain Parameters of Equation, Mesh and Solver form user defined values

prob = 1;       %  1: Conv-Diff, 2: Poisson
   N = 30;      %  number of intervals of the grid in x-dir (must be even)
   M = 40;      %  number of intervals of the grid in y-dir (must be even)
epsi = 1e-4;    %  scalar diffusion coefficient (perturbation parameter)

[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);

%% Generate mesh

[MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );

%% Construct Linear System

%  Construct matrix of linear system using Kroenecker products
[SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);

%  Obtain RHS of linear system by applying contributions of B.C.'s
[SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);

%% Solve system by direct inversion and GMRES method

[SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );

%% Solve system by the mult. Schwarz method with T = Q2Q1 and T = Q1Q2

SOLVER.order = 1;

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

%% Plot Results

%  Algebraic error generated by the Schwarz method
block_Schwarz_ploterr(EQ, SOLUTION);

%  Solutions obtained by the different approaches
[max_err] = block_Schwarz_plotsol( MESH, SOLUTION);

problem_size = size(SOLVER.A,1)
max_err
