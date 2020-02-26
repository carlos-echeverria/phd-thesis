function [ SOLUTIONparams ] = block_Schwarz_direct_solve(SOLVERparams, MESHparams, EQparams)
%BLOCK_SCHWARZ_GMRES uses the GMRES method implemented by MATLAB to compute
%the solution to the system Ax=b, where A comes from the Shishkin
%discretization of the 2D convection diffusion equation and b is the
%discretized right hand side of the following BVP:
%
%  -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma.
%
%   When there exists an analytic solution to the above problem, it is also
%   evaluaed on the grid points of the mesh.
%
%   function call: 
%
%     [SOLUTIONparams] = block_Schwarz_direct_solve(SOLVERparams)
%
%   input:
%
%    SOLVERparams: data structure with solver parameters including matrix 
%                  'A' and rhs 'b'. 
%
%   output:
%
%  SOLUTIONparams: structure array with the following solver parameters:
%                  x_g: solution of linear system obtained by GMRES. 
%                res_g: vector with the residulas at each iteration.
%                  x_d: vector with the solution of linear system obtained 
%                       by backslash operator.
%                 u_ex: matrix with the analytic solution evaluated at the 
%                       mesh points (including boundary points)
%                 x_ex: vector with the analytic solution evaluated at the
%                       internal mesh points (excluding boundary). 
%
%
% Written by Carlos Echeverria on August 8, 2018.
% Last Edited by C. E. on September 24, 2019.


  tol = SOLVERparams.tol;      
maxit = SOLVERparams.maxit;  
   x0 = SOLVERparams.x0; 
    A = SOLVERparams.A;   
    b = SOLVERparams.b;

    u = EQparams.u;
    
    N = MESHparams.N;
    M = MESHparams.M;
   xi = MESHparams.xi;
   yi = MESHparams.yi;
  
  tic; 
  x_d = A\b;
  time_direct=toc;
  x_g = x0;

tic;
[x_g,~,~,~,res_g] = gmres(A,b,[],tol,maxit,[],[],x_g);
time_gmres=toc;


   SOLUTIONparams.x_d = x_d;
   SOLUTIONparams.x_g = x_g;
 SOLUTIONparams.res_g = res_g;
 SOLUTIONparams.time_g = time_gmres;
 
 %  Evaluate exact solution on the grid and store it
     [X,Y] = meshgrid(xi,yi); 
      u_ex = u(X',Y');
   SOLUTIONparams.u_ex = u_ex;
   
      x_eX = u_ex(2:N,2:M);
      x_ex = reshape(x_eX,[(N-1)*(M-1),1]);
      
   SOLUTIONparams.x_ex = x_ex;
   SOLUTIONparams.time_d = time_direct;

end