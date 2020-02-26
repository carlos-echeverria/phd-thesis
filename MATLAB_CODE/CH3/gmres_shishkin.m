function [ ] = gmres_shishkin(epsilon, alpha, beta, N)
%GMRES_SHISHKIN performs numerical examples for the behavior of GMRES
% for the Shishkin mesh discretized 1D convection diffusion problem:
%
% -epslion u'' + alpha u' + beta u = 1, in (0,1), u(0) = u(1) = 0
%
%   For each set of parameters epsilon, alpha, beta, and N, the function
%   produces a Figure that compares the residual reduction history of the
%   GMRES method when the matrix is obtained via different discretization
%   schemes. The upwind as well as the central diference scheme are
%   compared together with their scaled versions.
%
%   function call:
%
%            [ ] = gmres_shishkin( epsilon, alpha, beta, N )
%
%    input:
%
%        epsilon: scalar diffusion coefficient (perturbation parameter)
%          alpha: scalar convection coefficient (wind)
%           beta: scalar reaction coefficient
%              N: number of intervals to subdivide the domain [0,1]
%
%   output:
%
%        Figure
%
%  subordinate functions:
%          
%        Shishkin_get_A.m
%
% Written by Carlos Echeverria on October 10, 2019.
% Last edited by C.E. on February 20, 2020.

%% Solver parameters

tol   = 1e-14;      % GMRES tolerance
maxit = N-1;        % maximum number of GMRES iterations
x0 = zeros(N-1,1);  % initial approximation

%% Experiments for upwind unscaled

disc = 1; scaled = 0;
[ A,~ ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = ones(N-1,1);
[~,~,~,~,res1] = gmres(A,b,[],tol,maxit,[],[],x0);

%% Experiments for upwind upwind scaled
disc = 1; scaled = 1;
[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = D*ones(N-1,1);
[~,~,~,~,res2] = gmres(A,b,[],tol,maxit,[],[],x0);

%% Experiments for central unscaled

disc = 2; scaled = 0;
[ A, ~ ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = ones(N-1,1);
[~,~,~,~,res3] = gmres(A,b,[],tol,maxit,[],[],x0);

%% Experiments for central scaled

disc = 2; scaled = 1;
[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = D*ones(N-1,1);
[~,~,~,~,res4] = gmres(A,b,[],tol,maxit,[],[],x0);

%% Produce figure

figure(1), clf;
purple=[.49 .18 .55];
blue=[.301 .745 .933];
semilogy(0:length(res1)-1,res1./res1(1),'k-', 'LineWidth', 2);
hold on;
% semilogy(0:length(res2)-1,res2./res2(1),'k:','LineWidth', 2);
semilogy(0:length(res3)-1,res3./res3(1),'r-.', 'LineWidth', 2);
% semilogy(0:length(res4)-1,res4./res4(1),'r--', 'LineWidth', 2);
%axis([0 12 1e-20 1.5]);
%set(gca,'XTick',[0:15]')
set(gca,'FontSize',18);
%hold off;
xlabel('k','FontSize',17);
ylabel('Residual norms','FontSize',17);
%title('GMRES','FontSize',16);
leg=legend('upwind','central','Location','SouthWest');
% leg=legend('upwind unscaled','upwind scaled','central unscaled','central scaled','Location','SouthWest');
set(leg,'FontSize', 17);
