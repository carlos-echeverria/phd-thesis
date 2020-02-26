function [] = experiments_gmres_precond_inexact_upwind(prob, N, M, epsi, local_tol)
%EXPERIMENTS_GMRES_PRECOND_INEXACT_UPWIND performs numerical experiments 
%   for solving the system: Ax = b, where the matrix A is obtained from the 
%   discretization of the following 2D conv-diff problem posed on a 
%   Shishkin mesh with one boundary layer near the outflow boundary:
%
%  -eps*Delta(u) + (0,1).Nabla(u) + u = f in Omega,  u = g on Gamma,
%
%   The system is solved via the GMRES method preconditioned with the 
%   *inexact* multiplicative Schwarz method. For each set of parameters:
%   epsi, N and M, the function produces a Figure that compares the 
%   residual reduction history of the iterates produced by the solution 
%   approach and compares it to the case of exact preconditioning. 
%
%   function call:
%
%       [ ] = experiments_gmres_precond_upwind(prob,N,M,epsi, local_tol)
%
%    input:
%
%      prob: type of problem. 1 for conv-diff, 2 for Poisson.
%         N: number of intervals of the Shishkin mesh in the x-direction.
%         M: number of intervals of the Shishkin mesh in the y-direction.
%      epsi: scalar diffusion coefficient (perturbation parameter).
% local_tol: desired accuracy of solution of local problems.
%
%   output:
%
%      * Figure with residual history of the preconditioned GMRES method
%        with exact and inexact preconditioning.
%      * Outputs values of number of inner and outer iterations, solution
%        time and error for different solution approaches.
%
%
% Written by Carlos Echeverria on October 10, 2019.
% Edited  by C.E. on February 25, 2020.

%% Define problem and obtain preconditioned system  (T and c) with exact localsolves 

[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);
            [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );
          [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);
          [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
        [SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );
        
        SOLVER.local_tol = local_tol;
        SOLVER.inexact = 0;
        [SOLUTION] = block_compute_T_c_inexact(SOLVER, MESH, SOLUTION);

   P1 = SOLUTION.P;
   T1 = SOLUTION.T;
   c1 = SOLUTION.c;

%    format longe
%        siz = size(T);
%        ran = rank(P);
%        sin = svd(P);
% %    pause

%% Solve preconditioned system with gmres

   tol   = 1e-7;      % GMRES tolerance
   maxit = ((N-1)*(M-1))-1;      % maximum number of GMRES iterations
   x0 = zeros(((N-1)*(M-1)),1);  % initial approximation

tic
[x1,~,~,~,res1] = gmres(P1,c1,[],tol,maxit,[],[],x0);
time_ex=toc;

%% Define problem and obtain preconditioned system  (T and c) with inexact localsolves 

[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);
            [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );
          [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);
          [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
        [SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );
        
        SOLVER.local_tol=local_tol;
        SOLVER.inexact=1;
        SOLUTION.resi = 0;
        [SOLUTION] = block_compute_T_c_inexact(SOLVER, MESH, SOLUTION);

      P2 = SOLUTION.P;
      T2 = SOLUTION.T;
      c2 = SOLUTION.c;
      
%% Solve preconditioned system with gmres

   tol   = 1e-7;                % GMRES tolerance
   maxit = ((N-1)*(M-1))-1;      % maximum number of GMRES iterations
   x0 = zeros(((N-1)*(M-1)),1);  % initial approximation

tic
[x2,~,~,~,res2] = gmres(P2,c2,[],tol,maxit,[],[],x0);
time_inex=toc;

%% Explicitly obtain T to calculate singular values (numerical rank)

%         I  = eye(size(A));
%        
%         for i=1:size(A,2)
%             
%             ei = I(:,i);
%             [~, c2,PP1(:,i),PP2(:,i), iter_h, res_h, iter_H, res_H] = compute_T(R1, R2, A, b, ei, order, local_tol);
% 
% %             PP1(:,i) = P1(ei);
% %             PP2(:,i) = P2(ei);
% 
%         end

%% plot convergence results
figure(1), clf
semilogy(0:length(res1)-1,res1./res1(1),'k-','LineWidth',2); hold on;
semilogy(0:length(res2)-1,res2./res2(1),'k:','LineWidth',2);
axis([0 10 1e-8 1e0]);
set(gca,'XTick',[0:1:10,12:2:26]')
set(gca,'FontSize',18)
xlabel('k','FontSize',17);
ylabel('Residual norms','FontSize',17);
title(sprintf('local solve tol= %5.1e', local_tol),'FontSize',17);
leg=legend('exact local solves','inexact local solves','Location','NorthEast');
set(leg,'FontSize',17)


%% output iterations and time to the console
fprintf('epsilon:%5.1e, tolerance:%5.1e, iter.exact(outer):%d, iter.inexact(outer):%d \n', epsi, local_tol, length(res1)-1, length(res2)-1);
fprintf('time (direct):%f sec., time(gmres):%f sec., time(exact schwarz):%f sec., time(inexact schwarz):%f sec. \n', SOLUTION.time_d, SOLUTION.time_g, time_ex, time_inex);


%% plot solution results
   N = MESH.N;
   M = MESH.M;
  xi = MESH.xi;
  yi = MESH.yi;
 x_d = SOLUTION.x_d;
 x_g = SOLUTION.x_g;
 u_ex = SOLUTION.u_ex;

fprintf('2-norm solution error. (dir.-gmres.)/(dir.) :%5.1e, (dir.-ex.)/(dir.) :%5.1e, (dir.-inex.)/(dir.) :%5.1e, \n', norm(x_d-x_g)/norm(x_d), norm(x_d-x1)/norm(x_d),norm(x_d-x2)/norm(x_d));

    [X,Y] = meshgrid(xi,yi); 

%  Embed solutions with boundary conditions:
  X_d  = u_ex; X_d(2:N, 2:M)    = reshape(x_d,  [N-1,M-1]); % direct solution (backslash)
  X_ex = u_ex; X_ex(2:N, 2:M)   = reshape(x1, [N-1,M-1]);   % gmres for preconditioned sys (exact localsolves)
X_inex = u_ex; X_inex(2:N, 2:M) = reshape(x2, [N-1,M-1]);   % gmres for preconditioned sys (inexact localsolves)

colormap('default')
figure(202), 
% set(gcf, 'Position',  [0, 0, 400, 700]);
    subplot(1,3,[1]), surf(X',Y',X_d), title('backslash solution','FontSize',17)
    axis([0 1 0 1 -1 1.0]); % view(90,0);
    set(gca,'FontSize',18)
    xlabel('x','FontSize',17);
    ylabel('y','FontSize',17);
    zlabel('u(x,y)','FontSize',17);
figure(202), 
    subplot(1,3,[2]), surf(X',Y',X_ex), title('precond. GMRES (exact local solves)','FontSize',17)
    axis([0 1 0 1 -1 1.0]); %view(90,0);
    set(gca,'FontSize',18)
    xlabel('x','FontSize',17);
    ylabel('y','FontSize',17);
    zlabel('u(x,y)','FontSize',17);
figure(202), 
    subplot(1,3,[3]), surf(X',Y',X_inex), title('precond. GMRES (inexact local solves)','FontSize',17)
     axis([0 1 0 1 -1 1.0]); %view(90,0);
     set(gca,'FontSize',18)
    xlabel('x','FontSize',17);
    ylabel('y','FontSize',17);
    zlabel('u(x,y)','FontSize',17);
     
     