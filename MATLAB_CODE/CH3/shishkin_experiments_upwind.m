function [] = shishkin_experiments_upwind(epsilon,alpha,beta,N,disc,scaled)
%   Perfroms the numerical experiments for solving the system:
%   Ax = b, where the matrix A is obtained from the upwind finit difference 
%   discretization of the following 1D conv-diff problem posed on a 
%   Shishkin mesh with one transition point:
%
%   -epsilon * u'' + alpha * u' + beta * u = f in (0,1),  u(0) = u(1) = 0.
%
%   The system is solved using the self implemented Multiplicative Schwarz 
%   method
%
%   function call:
%
% [ A,D ] = shishkin_experiments_upwind(epsilon,alpha,beta,N,disc,scaled)
%
%   input: 
%         epsilon: scalar diffusion  coefficient.  
%           alpha: scalar convection coefficient.
%            beta: scalar reaction   coefficient.
%               N: the number of intervals conformng the grid.
%            disc: specifies the type of discretization.
%                  can be 1 for upwind or 2 for central differencing.
%          scaled: can be 1 for a scaled matrix or 0 for unscaled.
%
%   output:
%               Figure: error curves with bounds
%              Console: numerical value of convergence factor and bound
%
%  subordinate functions:
%          
%        Shishkin_get_A.m
%        Multiplicative_Schwarz.m
%
% Written by Carlos Echeverria on September 10, 2016.
% Edited  by Joerg Liesen, September 26, 2016
% Edited  by C.E. on February 20, 2020.


tau = 2*epsilon*log(N)/alpha;
n = N/2; 
H = (1-tau)/n; 
h = tau/n;
m = (N/2) - 1;

conv_factor = epsilon / ( epsilon + alpha*H );

%% Other parameters
tol   = 1e-14;   % tolerance for the mSm method
maxit = 30;      % maximum number of iterations for the mult. Schwarz method

[ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = D*ones(N-1,1);
x = A\b;
x0 = zeros(N-1,1);

% order 1: T=(I-P2)(I-P1)
[ ~, ~, err_s1, ~, T] = Multiplicative_Schwarz( A,b,maxit,tol,1,x0,x,0);

rho = abs(T(n+1,n+1));
fprintf('N = %d, epsilon = %5.1e, scheme =  upwind, rho = %5.1e, bound = %5.1e \n',N, epsilon,rho,conv_factor);

bound1(1) = 1;
for j = 1:length(err_s1)-1
    bound1(j+1) = conv_factor * bound1(j);
end

% order 2: T=(I-P1)(I-P2)
[ ~, ~, err_s2] = Multiplicative_Schwarz( A,b,maxit,tol,2,x0,x,0);
bound2(1) = 1; bound2(2) = 1;
for j = 2:length(err_s2)-1
    bound2(j+1) = conv_factor * bound2(j);
end

figure(1), clf
semilogy(0:length(err_s1)-1,err_s1/err_s1(1),'k-.', 'LineWidth', 2); hold on;
semilogy(0:length(bound1)-1,bound1,'rs', 'MarkerSize', 10);
semilogy(0:length(err_s2)-1,err_s2/err_s2(1),'k-', 'LineWidth', 2);
semilogy(0:length(bound2)-1,bound2,'r^', 'MarkerSize', 10);
axis([0 30 1e-16 1e0]);
set(gca,'XTick',[0:3:30]')
set(gca,'FontSize',18)
xlabel('k','FontSize',17);
ylabel('Error norms and bounds','FontSize',17);
title('Upwind','FontSize',17);
leg=legend('T=Q_2Q_1','error bound','T=Q_1Q_2','error bound','Location','Best');
set(leg,'FontSize',17);
