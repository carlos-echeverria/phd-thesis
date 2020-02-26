function [] = block_Schwarz_plotsol_inexact(MESHparams, SOLUTIONparams, x_ex, x_inex)
%BLOCK_SCHWARZ_PLOTSOL plots the solution to the boundary value problem:
%
% -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma,
%
%  on the domain Omega = [0,1]x[0,1] discretized by a hybrid Shishkin 
%  mesh with N intervals on the x-direction and M intervals on the 
%  y-direction. The solutions presented are calculated by a direct solution
%  method (MATLAB's backslash) and an iterative method (Multiplicatie Schwarz). 
%  The analytic solution is also presented when available.
%
%   function call:
%
%        [ ] = block_Schwarz_plotsol(MESHparams,SOLUTIONparams)
%
%   input: 
%
%      MESHparams: data structure the following mesh parameters:
%              N: even integer with the number of intervals in the x-dir. 
%              M: even integer with the number of intervals in the y-dir.
%             xi: vector with the coordinates of the points in the x-dir.
%             yi: vector with the coordinates of the points in the y-dir.
%
%  SOLUTIONparams: data structure with the following solution parameters:
%            x_d: (N-1)*(M-1) vector with solution obtained with 
%                  MATLAB's backsalsh operator. 
%            x_s1: (N-1)*(M-1) vector with solution obtained with 
%                   with the multiplicative Schwarz method.
%            u_ex: (N-1)*(M-1) vector with solution obtained by
%                   evaluating the analytic solution on the grid points.
%
%   output:
%
%             figure plot of the different solutions.
%             figure plot of the global error distribution.
%
%
% Written by Carlos Echeverria on September 24, 2016.
% Last Edited by C. E. on September 24, 2019.


   N = MESHparams.N;
   M = MESHparams.M;
  xi = MESHparams.xi;
  yi = MESHparams.yi;
 x_d = SOLUTIONparams.x_d;
%  x_s = SOLUTIONparams.x_s;
% x_s2 = SOLUTIONparams.x_s;
% x_g = SOLUTIONparams.x_g;
% x_ex = SOLUTIONparams.x_ex;
u_ex = SOLUTIONparams.u_ex;

 
%% Plot Solutions on the grid

    [X,Y] = meshgrid(xi,yi); 
  % [X,Y] = meshgrid(1:length(xi),1:length(yi));

%  Embed solutions with boundary conditions:
  X_d  = u_ex; X_d(2:N, 2:M)  = reshape(x_d,  [N-1,M-1]); % direct solution (backslash)
  X_ex = u_ex; X_ex(2:N, 2:M) = reshape(x_ex, [N-1,M-1]); % Schwarz solution (T=Q2Q1)
X_inex = u_ex; X_inex(2:N, 2:M) = reshape(x_inex, [N-1,M-1]); % Schwarz solution (T=Q2Q1)

%  X_g = u_ex;  X_g(2:N, 2:N)  = reshape(x_g, [N-1,N-1]); % GMRES solution


colormap('default')
figure(202), 
% set(gcf, 'Position',  [0, 0, 400, 700]);
    subplot(3,2,[1, 2]), surf(X',Y',X_d), title('Analytic solution')
    axis([0 1 0 1 -1 2.5]); % view(90,0);
figure(202), 
    subplot(3,2,[3, 4]), surf(X',Y',X_ex), title('exact local solves (backslash)')
    axis([0 1 0 1 -1 2.5]); %view(90,0);
figure(202), 
    subplot(3,2,[5, 6]),surf(X',Y',X_inex), title('inexact local solves (gmres)')
    axis([0 1 0 1 -1 2.5]); %view(90,0);
    

end