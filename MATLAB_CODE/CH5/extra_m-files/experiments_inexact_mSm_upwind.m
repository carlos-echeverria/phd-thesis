function [] = experiments_inexact_mSm_upwind(prob, N, M, epsi, local_tols)

[EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);
            [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );
          [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);
          [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
        [SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ );
  SOLVER.order = 1; % order 1: T=(I-P2)(I-P1)
        [SOLUTION] = block_Multiplicative_Schwarz(SOLVER, MESH, SOLUTION);

% store computed solution, error, and residual of exact method.
err_s = SOLUTION.err_s_inf;
  x_s = SOLUTION.x_s;
res_s = SOLUTION.res_s;

figure(101), clf;
    semilogy(0:length(err_s)-1, err_s/err_s(1), 'k-' , 'LineWidth', 4); hold on


for j = 1:length(local_tols)

    SOLVER.local_tol = local_tols(j);
    [SOLUTION] = block_Multiplicative_Schwarz_inexact(SOLVER, MESH, SOLUTION);

    % stores computed solution, error and residual of inexact methods
     err_s2(:,j) = SOLUTION.err_s_in;
       x_s2(:,j) = SOLUTION.x_s_in;
     res_s2(:,j) = SOLUTION.res_s_in;

     % plot results
     figure(101),
     linS = {':','-.','--',':','-.','--',':','-.','--',':','-.','--'};
     semilogy(0:length(err_s2(:,j))-1, err_s2(:,j)/err_s2(1,j), 'LineStyle' , linS{j}, 'LineWidth', 4); hold on

end

set(gca,'FontSize',18);
xlabel('Iteration step: k','FontSize',17);
ylabel('Error norms: ||e^{(k)}||_{\infty}','FontSize',17);
title(['mSm exact/inexact, epsilon = ' num2str(epsi,'%1.0e')],'FontSize',17);
axis([0 10 1e-15 1e1]);
set(gca,'XTick',[0:1:10]')
leg = legend('exact', ['inexact' num2str(local_tols(1),'%1.0e')], ['inexact' num2str(local_tols(2),'%1.0e')],['inexact' num2str(local_tols(3),'%1.0e')],['inexact' num2str(local_tols(4),'%1.0e')],['inexact' num2str(local_tols(5),'%1.0e')],['inexact' num2str(local_tols(6),'%1.0e')],'Location', 'SouthEast');
% leg = legend('exact', ['inexact' num2str(local_tols(1),'%1.0e')], ['inexact' num2str(local_tols(2),'%1.0e')],['inexact' num2str(local_tols(3))],['inexact' num2str(local_tols(4))],'Location', 'SouthEast');
set(leg,'FontSize',17);


% code below can be uncomented to compare solutions and algebraic error:

% Plot solutions obtained by the different approaches on the grid
%    N = MESH.N;
%    M = MESH.M;
%   xi = MESH.xi;
%   yi = MESH.yi;
%  x_d = SOLUTION.x_d;
% u_ex = SOLUTION.u_ex;
%
% [X,Y] = meshgrid(xi,yi);
%
% %  Embed solutions with boundary conditions:
% X_d  = u_ex;  X_d(2:N, 2:M) = reshape(x_d, [N-1,M-1]); % backslash sol
% X_s1 = u_ex; X_s1(2:N, 2:M) = reshape(x_s, [N-1,M-1]); % exact schwarz sol
% X_s2 = u_ex; X_s2(2:N, 2:M) = reshape(x_s2(:,1), [N-1,M-1]); % Schwarz sol
% %  X_g = u_ex;  X_g(2:N, 2:N)  = reshape(x_g, [N-1,N-1]); % GMRES solution
%
%  colormap('default')
%  figure(202), set(gcf, 'Position',  [0, 0, 400, 700]);
%      subplot(3,2,[1, 2]), surf(X',Y',u_ex), title('Direct solution (backslash)')
%      axis([0 1 0 1 -1 2.5]); % view(90,0);
%  figure(202),
%      subplot(3,2,[3, 4]), surf(X',Y',X_d), title('Exact Schwarz solution (T=Q_2Q_1)')
%      axis([0 1 0 1 -1 2.5]); %view(90,0);
%  figure(202),
%      title_fig = sprintf('Inexact Schwarz solution local tol.=%0.0e',local_tols(1));
%      subplot(3,2,[5, 6]),surf(X',Y',X_s1), title(title_fig)
%      axis([0 1 0 1 -1 2.5]); %view(90,0);
%
%
% %% Plot Error on the grid
%               error  = abs(x_s2(:,1)-x_s);
%              max_err = max(error);
%               Error  = zeros(N+1,M+1);
%     Error(2:N, 2:M)  = reshape(error,[N-1,M-1]);
%
% figure(303), set(gcf, 'Position',  [500, 20, 700, 300]);
% title_f = sprintf('Inexact Schwarz solution (local tol.=%5.0e)',local_tols(1));
%    subplot(2,2,[1,3]),surf(X',Y',X_s1),title(title_f)
% figure(303),
%    subplot(2,2,[2,4]),surf(X',Y',Error),title('Error distribution (w.r.t exact Schwarz sol.)')
%

end
