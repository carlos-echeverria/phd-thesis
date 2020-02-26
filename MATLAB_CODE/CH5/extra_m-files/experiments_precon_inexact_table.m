function [res1,res2]= experiments_precon_inexact_table(prob,N,M,epsi,local_tol)
%GMRES_SHISHKIN_PRECON_INEXACT_TABLE performs numerical examples for the
% behavior
% of the preconditioned GMRES methd with the multiplicative Schwarz method
% as a preconditioner for the Shishkin mesh discretized 1D convection
% diffusion problem:
%
%     -epslion u'' + alpha u' + beta u = 1, in (0,1), u(0) = u(1) = 0
%
%   For each set of parameters epsilon, alpha, beta, and N,
%   the function produces a Figure that compares the residual reduction
%   history of the GMRES method preconditioned with the multilpicative
%   Schwarz method when the matrix is obtained via different discretization
%   schemes. The upwind as well as the central diference scheme are
%   compared together with their scaled versions.
%
%   function call:
%
% [ ] = gmres_shishkin_precon_inexact(epsilon, alpha, beta, N,order,local_tol)
%
%    input:
%
%        epsilon: scalar diffusion coefficient (perturbation parameter)
%          alpha: scalar convection coefficient (wind)
%           beta: scalar reaction coefficient (wind)
%              N: number of intervals to subdivide the domain [0,1]
%          order:
%      local_tol:
%
%   output:
%
%        Figure
%
% Written by Joerg Liesen on September 14, 2016
% Modified for the preconditioner by Daniel B. Szyld on December 30, 2016
% modified for the thesis by Carlos Echeverria on October 16, 2019

%% Other parameters
tol   = 1e-14;      % GMRES tolerance
maxit = N-1;        % maximum number of GMRES iterations
x0 = zeros(N-1,1);  % initial approximation

%% upwind unscaled exact
%% upwind unscaled exact
[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);
            [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );
          [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);
          [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
        [SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );
        SOLVER.local_tol=local_tol;
        SOLVER.inexact=0;
        [SOLUTION] = block_compute_T_c_inexact(SOLVER, MESH, SOLUTION);

   P = SOLUTION.P;
   c = SOLUTION.c;

%% Other parameters
   tol   = 1e-14;      % GMRES tolerance
   maxit = ((N-1)*(M-1))-1;      % maximum number of GMRES iterations
   x0 = zeros(((N-1)*(M-1)),1);  % initial approximation
[~,~,~,~,res1] = gmres(P,c,[],tol,maxit,[],[],x0);

%% upwind unscaled inexact
[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);
            [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );
          [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);
          [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
        [SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );
        SOLVER.local_tol=local_tol;
        SOLVER.inexact=1;
        [SOLUTION] = block_compute_T_c_inexact(SOLVER, MESH, SOLUTION);

   P = SOLUTION.P;
   c = SOLUTION.c;

   %% Other parameters
   tol   = 1e-14;      % GMRES tolerance
   maxit = ((N-1)*(M-1))-1;      % maximum number of GMRES iterations
   x0 = zeros(((N-1)*(M-1)),1);  % initial approximation
[~,~,~,~,res2] = gmres(P,c,[],tol,maxit,[],[],x0);



%% plot results
figure(1), clf
semilogy(0:length(res1)-1,res1./res1(1),'k-','LineWidth',2); hold on;
semilogy(0:length(res2)-1,res2./res2(1),'k-.','LineWidth',2);
% semilogy(0:length(res3)-1,res3./res3(1),'r-','LineWidth',4);
% semilogy(0:length(res4)-1,res4./res4(1),'r-.','LineWidth',4);
axis([0 10 1e-16 1.5]);
set(gca,'XTick',[0:1:10]')
set(gca,'FontSize',18)
%hold off;
xlabel('k','FontSize',17);
ylabel('Residual norms','FontSize',17);
title(sprintf('inexact localsolves tol= %5.1e', local_tol),'FontSize',17);
%legend('upwind unscaled','upwind scaled','central unscaled','central scaled','Location','SouthWest');
leg=legend('upwind unscaled exact','upwind unscaled  inexact','central unscaled exact','central unscaled inexact','Location','SouthEast');
set(leg,'FontSize',17)
