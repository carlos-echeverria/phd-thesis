function [ ] = gmres_shishkin_precond(epsilon, alpha, beta, N)
%GMRES_SHISHKIN_PRECOND performs numerical examples for the behavior of
% GMRES for the Shishkin mesh discretized 1D convection diffusion problem:
%
% -epslion u'' + alpha u' + beta u = 1, in (0,1), u(0) = u(1) = 0
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
%          [ ] = gmres_shishkin_precond(epsilon, alpha, beta, N)
%
%    input:
%
%        epsilon: scalar diffusion coefficient (perturbation parameter)
%          alpha: scalar convection coefficient (wind)
%           beta: scalar reaction coefficient (wind)
%              N: number of intervals to subdivide the domain [0,1]
%
%  subordinate functions:
%          
%        Shishkin_get_A.m
%        Compute_T_c.m
%
%   output:
%
%        Figure
%
% Written by Carlos Echeverria on September 10, 2016.
% Written by Joerg Liesen on September 14, 2016
% Modified for the preconditioner by Daniel B. Szyld on December 30, 2016
% Edited by Carlos Echeverria on October 16, 2019
% Edited  by C.E. on February 20, 2020.

%% Other parameters
tol   = 1e-14;      % GMRES tolerance
maxit = N-1;        % maximum number of GMRES iterations
x0 = zeros(N-1,1);  % initial approximation

%% upwind unscaled
disc = 1; scaled = 0;
[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = ones(N-1,1);
[ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c( A,b,disc);
[~,~,~,~,res1] = gmres(B,c,[],tol,maxit,[],[],x0);

%% upwind scaled
disc = 1; scaled = 1;
[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = D*ones(N-1,1);
[ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c( A,b,disc);
[~,~,~,~,res2] = gmres(B,c,[],tol,maxit,[],[],x0);

%% central unscaled
disc = 2; scaled = 0;
[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = ones(N-1,1);
[ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c( A,b,disc);
[~,~,~,~,res3] = gmres(B,c,[],tol,maxit,[],[],x0);

%% central scaled
disc = 2; scaled = 1;
[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = D*ones(N-1,1);
[ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c( A,b,disc);
[~,~,~,~,res4] = gmres(B,c,[],tol,maxit,[],[],x0);

%% Produce Figure
figure(1), clf
semilogy(0:length(res1)-1,res1./res1(1),'k-', 'LineWidth',2); hold on;
semilogy(0:length(res2)-1,res2./res2(1),'k:', 'LineWidth',2);
semilogy(0:length(res3)-1,res3./res3(1),'r-.', 'LineWidth',2);
semilogy(0:length(res4)-1,res4./res4(1),'r--', 'LineWidth',2);
axis([0 10 1e-16 1e0]);
set(gca,'XTick',[0:1:10]')
set(gca,'FontSize',18)
%hold off;
title('exact local solves','FontSize',17)
xlabel('k','FontSize',17);
ylabel('Residual norms','FontSize',17);
%title('GMRES','FontSize',16);
%legend('upwind unscaled','upwind scaled','central unscaled','central scaled','Location','SouthWest');
leg=legend('upwind unscaled','upwind scaled','central unscaled','central scaled','Location','NorthEast');
set(leg,'FontSize',17)
