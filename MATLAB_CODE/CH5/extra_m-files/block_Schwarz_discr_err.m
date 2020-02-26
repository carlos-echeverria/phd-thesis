function [globalerr, globalerr2, max_err_vec ] = block_Schwarz_discr_err( prob,epsi )
%BLOCK_SCHWARZ_PLOTSOL plots the solution to the boundary value 
% problem:
%
% -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma,
%
%  on the domain Omega = [0,1]x[0,1] discretized by a hybrid Shishkin 
%  mesh with N intervals on each coordinate direction. The solutions
%  presented are calculated by a direct solution method (MATLAB's backslash)
%  and an iterative method (Multiplicatie Schwarz). The analytic solution
%  is also presented when availabel.
%
%   function call:
%
% [SOLUTIONparams] = block_Schwarz_plotsol(EQparams,MESHparams,SOLUTIONparams)
%
%   input: 
%             N: number of intervals in each oordinate direction.
%            xi: column vector with x-coordinates of the mesh.
%            yi: column vector with y-coordinates of the mesh.
%             u: matrix with analtic solution in each entry.
%           x_d: direct or "exact" solution vector x_d=A\b.
%           x_s: solution vector by the multiplicative Schwarz method. 
%                (obtained with Multiplicative_Schwarz_block_bounds.m)
%           x_g: solution vector  obtained by GMRES.     
%             f: discrete right hand side of the PDE with boundary cond.
%
%   output:
%             figure plot of the different solutions.
%
%
% Written by Carlos Echeverria on January 15, 2019.
% Last Edited by C. E. on September 24, 2019.
    

     npts = 2.^(3:1:6);

globalerr = zeros(length(npts),1);
globalerr2 = zeros(length(npts),1);
max_err_vec = zeros(length(npts),length(npts));

for ii = 1:length(npts)
    for jj = 1:length(npts)

    
     M = npts(ii);
     N = npts(jj); % warning: keep points in one direction relatively small
  
  
    [EQ, MESH, SOLVER] = block_Schwarz_getparams(prob, N, M, epsi);
                [MESH] = block_Schwarz_2D1L_mesh( MESH, EQ );
              [SOLVER] = block_Schwarz_2D1L_get_b(SOLVER, MESH, EQ);
              [SOLVER] = block_Schwarz_2D1L_Kron_get_A(SOLVER,MESH,EQ);
            [SOLUTION] = block_Schwarz_direct_solve(SOLVER, MESH, EQ ); 
          SOLVER.order = 1; % choose order of solution of subproblems 1:Q2Q1  2:Q1Q2
            [SOLUTION] = block_Multiplicative_Schwarz(SOLVER, MESH, SOLUTION);
            [max_err]  = block_Schwarz_plotsol( MESH, SOLUTION);
                    
                  x_d  = SOLUTION.x_d;  % direct solution (backslash)
                   x_s = SOLUTION.x_s;  % Schwarz solution (T=QiQj)
                  x_ex = SOLUTION.x_ex; % exact solution (analytic)

       max_err_vec(ii,jj) = max_err;                  
        globalerr(jj)  = norm(x_s-x_ex, Inf)/norm(x_ex,Inf);
        globalerr2(jj) = norm(x_d-x_ex, Inf)/norm(x_d,Inf);
        
        clear EQ MESH SOLVER SOLUTION
    
    end
end
    
figure(102), loglog(npts, globalerr, 'MarkerSize', 10), grid on,
title('Global Error','FontSize',16)
set(gca,'FontSize',18)
xlabel('Number of grid points in refined direction','FontSize',16);
ylabel('Relative Error Norm (||x_s-x_d||/||x_d||)' ,'FontSize',16);

figure(103), loglog(npts, globalerr2, 'MarkerSize', 10), grid on,
title('Global Error2','FontSize',16)
set(gca,'FontSize',18)
xlabel('Number of grid points in refined direction','FontSize',16);
ylabel('Relative Error Norm (||x_s-x_d||/||x_d||)' ,'FontSize',16);

figure(104), plot(npts, max_err_vec, 'MarkerSize', 10), grid on,
title('MAX Error','FontSize',16)
set(gca,'FontSize',18)
xlabel('Number of grid points in refined direction','FontSize',16);
ylabel('Relative Error Norm (||x_s-x_d||/||x_d||)' ,'FontSize',16);

end