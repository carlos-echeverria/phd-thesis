function [] = gmres_shishkin_precon_inexact(epsilon,alpha,beta,N,order,local_tols)
%GMRES_SHISHKIN_PRECON_INEXACT performs numerical examples for the behavior
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
% Written by Carlos Echeverria on September 10, 2016.
% Written by Joerg Liesen on September 14, 2016
% Modified for the preconditioner by Daniel B. Szyld on December 30, 2016
% Edited by Carlos Echeverria on October 16, 2019
% Edited  by C.E. on February 20, 2020.

%% Other parameters
tol   = 1e-14;      % GMRES tolerance
maxit = N-1;        % maximum number of GMRES iterations
x0 = zeros(N-1,1);  % initial approximation

%% upwind exact
disc = 1; % upwind
scaled = 0; % unscaled
[ A,~] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
b = ones(N-1,1);
[ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c( A, b, order);
[x_ex,~,~,~,res1] = gmres(B,c,[],tol,maxit,[],[],x0);

 
figure(111), clf % plot residuals
    semilogy(0:length(res1)-1,res1./res1(1),'k-','LineWidth',2); hold on;

n = N/2;
tau = min(0.5, 2*eps*log(N)/alpha);
H = 2*(1-tau)/N; h = 2*tau/N;

for i=1:n
    mesh(i)=i*H;
end

for i=n+1:N-1
    mesh(i)=1-(N-i)*h;
end
    
figure(112), clf % plot solution 
%     plot(mesh,x_ex,'k-','LineWidth',2); hold on;
    plot(x_ex,'k-','LineWidth',2); hold on;


%% upwind inexact
inexact=1;
for j=1:length(local_tols)
    local_tol=local_tols(j);
    [ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c_inexact( A,b,order, inexact, local_tol);
             [x_inex,~,~,~,res_inex] = gmres(B,c,[],tol,maxit,[],[],x0);
    
    fprintf('epsilon:%5.1e, inexact tol:%5.1e, upwind (exact):%d, upwind (inexact):%d\n', epsilon, local_tol, length(res1)-1, length(res_inex)-1);

    figure(111),
     linS = {':','-.','--',':','-.','--',':','-.','--',':','-.','--'};
     semilogy(0:length(res_inex)-1,res_inex./res_inex(1),'LineStyle' , linS{j}, 'LineWidth', 2); hold on
    figure(112),
%      plot(mesh,x_inex,'LineStyle' , linS{j}, 'LineWidth', 2); hold on
     plot(x_inex,'LineStyle' , linS{j}, 'LineWidth', 2); hold on

end

figure(111),
axis([0 10 1e-16 1e0]);
set(gca,'XTick',[0:1:10]')
set(gca,'FontSize',18)
%hold off;
xlabel('k','FontSize',17);
ylabel('Residual norms','FontSize',17);
title(sprintf('inexact localsolves tol= %5.1e', local_tol),'FontSize',17);
leg=legend('upwind exact','upwind 1e-1','upwind 1e-2','Location','NorthEast');
set(leg,'FontSize',17)

     % %% central unscaled exact
% disc = 2; scaled = 0;
% [ A,~] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
% b = ones(N-1,1);
% [ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c( A,b,order);
% [~,~,~,~,res3] = gmres(B,c,[],tol,maxit,[],[],x0);
% 
% %% central  inexact
% disc = 2; scaled = 1;
% [ A,D ] = Shishkin_get_A( epsilon,alpha,beta,N,disc,scaled );
% b = D*ones(N-1,1);
% inexact=1;
% [ c, ~, B, ~, ~, ~, ~ ] = Compute_T_c_inexact( A,b,order, inexact, local_tol);
% [~,~,~,~,res4] = gmres(B,c,[],tol,maxit,[],[],x0);


% semilogy(0:length(res3)-1,res3./res3(1),'r--','LineWidth',2);
% semilogy(0:length(res4)-1,res4./res4(1),'r-.','LineWidth',2);
% axis([0 10 1e-16 1e0]);
% set(gca,'XTick',[0:1:10]')
% set(gca,'FontSize',18)
% %hold off;
% xlabel('k','FontSize',17);
% ylabel('Residual norms','FontSize',17);
% title(sprintf('inexact localsolves tol= %5.1e', local_tol),'FontSize',17);
% %legend('upwind exact','upwind scaled','central unscaled','central scaled','Location','SouthWest');
% leg=legend('upwind exact','upwind inexact','central exact','central inexact','Location','NorthEast');
% set(leg,'FontSize',17)

% fprintf('epsilon:%5.1e, inexact:%5.1e, upwind:%d, central:%d \n', epsilon, local_tol, length(res1)-1, length(res3)-1);

