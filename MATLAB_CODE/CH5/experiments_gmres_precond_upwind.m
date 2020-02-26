function [ ] = experiments_gmres_precond_upwind(prob,N,M,epsi )
%EXPERIMENTS_GMRES_PRECOND_UPWIND performs numerical experiments for 
%   solving the system: Ax = b, where the matrix A is obtained from the 
%   discretization of the following 2D conv-diff problem posed on a 
%   Shishkin mesh with one boundary layer near the outflow boundary:
%
%  -eps*Delta(u) + (0,1).Nabla(u) + u = f in Omega,  u = g on Gamma,
%
%   The system is solved via the GMRES method preconditioned with the mul-
%   tiplicative Schwarz method. For each set of parameters epsi, N and M,
%   the function produces a Figure that compares the residual reduction
%   history of the iterates produced by the solution approach. 
%
%   function call:
%
%          [ ] = experiments_gmres_precond_upwind(prob,N,M,epsi )
%
%    input:
%
%      prob: type of problem. 1 for conv-diff, 2 for Poisson.
%         N: number of intervals of the Shishkin mesh in the x-direction.
%         M: number of intervals of the Shishkin mesh in the y-direction.
%      epsi: scalar diffusion coefficient (perturbation parameter).
%
%   output:
%
%      * Figure with residual history of the preconditioned GMRES method
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

%% obtain iteration matrix T and correction vector c

[SOLUTION] = block_compute_T_c(SOLVER, MESH, SOLUTION);

   P = SOLUTION.P;
   c = SOLUTION.c;

%% set GMRES solution parameters
   tol   = 1e-14;      % GMRES tolerance
   maxit = ((N-1)*(M-1))-1;        % maximum number of GMRES iterations
   x0 = zeros((N-1)*(M-1),1);  % initial approximation
    [~,~,~,~,res1] = gmres(P,c,[],tol,maxit,[],[],x0);


%% Produce Figure
figure(111), clf,
semilogy(0:length(res1)-1,res1./res1(1),'k-', 'LineWidth',2); hold on;
% semilogy(0:length(res2)-1,res2./res2(1),'k:', 'LineWidth',2);
% semilogy(0:length(res3)-1,res3./res3(1),'r-.', 'LineWidth',2);
% semilogy(0:length(res4)-1,res4./res4(1),'r--', 'LineWidth',2);
axis([0 10 1e-16 1e0]);
set(gca,'XTick',[0:1:10]);
set(gca,'FontSize',18);
%hold off;
title('Preconditioned GMRES','FontSize',17)
xlabel('k','FontSize',17);
ylabel('Residual norms','FontSize',17);
%title('GMRES','FontSize',16);
%legend('upwind unscaled','upwind scaled','central unscaled','central scaled','Location','SouthWest');
leg=legend('upwind unscaled','upwind scaled','central unscaled','central scaled','Location','NorthEast');
set(leg,'FontSize',17)
