function [EQparams, MESHparams, SOLVERparams] = block_Schwarz_getparams(prob, N, M, epsi)
%BLOCK_SCHWARZ_GETPARAMS defines the numerical values of the parameters of 
% the governing equation for the following elliptic BVP:
%
%  -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma,
%
% where eps is the scalar diffusion coefficient, alpha = (wx, wy) is the 
% 'wind vector', and beta is the scalar reaction coefficient. 
%
% This function also defines the parameters of the matrix A resulting from
% the discretization of the previous equation. The parameters are
% implicitly given by defining the number of intervals and points of the 
% given mesh. Other parameters are also defined (see code).
%
% Lastly, the tolerance, maximum number of iterations, and initial value 
% of the iterative solver are also defined.
%
%    ALL OF THE INITIAL PARAMETER VALUES MUST BE HARDCODED BY THE USER!
%
%   function call: 
%
%[EQparams, MESHparams, SOLparams] = block_Schwarz_getparams(prob, N, epsi)
%
%    input: 
%
%         prob: scalar specifying problem - 1: ConvDiff, 2: Laplace.
%            N: Number of intervals in each direction of the mesh. It needs
%               to be even in order to have matrices with the desired
%               block structure.
%         epsi: scalar diffusion coefficient (perturbation parameter).
%
%   output:
%       
%      EQparams: data structure array with the following subfileds:
%               wx: component of convection coefficient in x direction.
%               wy: component of convection coefficient in y direction.
%             beta: scalar reaction coefficient.
%                f: function handle defining rhs of PDE.
%                g: function handle defining boundary function of BCs.
%                u: function handle defining analytic solution of PDE.
%                   (when available).
%              ph1: separable x-component of wx.
%              ps1: separable y-component of wx.
%              ph1: separable x-component of wy.
%              ps2: separable y-component of wy.
%
% The last four parameters are needed for the Kronecker formulation of the 
% matrix. In this case we we are restricted to using 'separable winds'
% that is, the components of alpha are separable functions in the space 
% variables, i.e., we need to define the convection vector as: 
%
%    alpha = (wx, wy) = ( phi1(x)*psi1(y), phi2(x)*psi2(y) ).
%
% For the case alpha = (0, 1) we have e.g., phi1=psi1=0, phi2=psi2=1.
%
%
%    MESHparams: data structure array with the following subfileds:
%                N: number of global intervals in x-direction.
%                M: number of global intervals in y-direction.
%           scaled: scaling parameter - 0:no, 1:yes.
%             disc: discretization parameter - 1:upwind, 2:central.
%
% The word 'internal' in the previous description refers to all the nodes 
% of the domain excluding the boundary nodes.
%              
%
% SOLVERparams: structure array with the following solver parameters:
%              tol: tolerance for the iterative method.
%            maxit: max number of iterations for the iterative method.
%               x0: initial approximation to the solution of the linsys.
%            order: order of application of the projections in the mult. 
%                   Schwarz iteration operator, 1: T = Q2Q1, 2: T = Q1Q2.
%
% Written by Carlos Echeverria on August 8, 2018.
% Edited  by C.E. on September 24, 2019.
% Edited  by C.E. on February  25, 2020.


EQparams.problem = prob; % 1: Conv-Diff, 2: Laplace

%% Definition of the EQparams data structure and its subfields
 
if prob==1;   % 1: Default parameters for the Conv-Diff problem
    
  EQparams.problem = prob; 
    EQparams.epsi  = epsi; 
      EQparams.wx  = 0;    
      EQparams.wy  = 1;    
    EQparams.beta  = 0;      
      EQparams.ph1 = 0; 
      EQparams.ps1 = 0; 
      EQparams.ph2 = 1; 
      EQparams.ps2 = 1;   
        EQparams.f = @(x,y) 0.*x+0.*y;
        EQparams.g = @(x,y) (2.*x-1).*(1-exp((y-1)/EQparams.epsi))/(1-exp(-1/EQparams.epsi));
        EQparams.u = @(x,y) (2.*x-1).*(1-exp((y-1)/EQparams.epsi))/(1-exp(-1/EQparams.epsi));


elseif prob==2; %2: Default parameters for the Laplace problem

  EQparams.problem = prob; 
    EQparams.epsi  = 1; 
      EQparams.wx  = 0;    
      EQparams.wy  = 0;    
    EQparams.beta  = 0;      
      EQparams.ph1 = 0; 
      EQparams.ps1 = 0; 
      EQparams.ph2 = 0; 
      EQparams.ps2 = 0; 
        EQparams.f = @(x,y) 4*((2*x.*y-1.25.*x-1.25.*y)./((x-1.25).^2+(y-1.25).^2));
        EQparams.g = @(x,y) -x.*y.*log((x-1.25).^2+(y-1.25).^2);  
        EQparams.u = @(x,y) -x.*y.*log((x-1.25).^2+(y-1.25).^2);
       
    % second set of parameters
%       EQparams.f = @(x,y) (2^2*pi^2+pi^2*4^2)*sin(2*pi*x).*sin(4*pi*y);
%       EQparams.g = @(x,y) 0.*x+0.*y;
%       EQparams.u = @(x,y) sin(2*pi*x).*sin(4*pi*y); 

end
     
   
%% Definition of the MESHparmas data structure and its subfields
  
     MESHparams.N  = N;        
     MESHparams.M  = M;       
 MESHparams.scaled = 1;
 MESHparams.disc   = 1;
 
 if mod(MESHparams.N, 2) ~= 0; 
     error('The number of intervals in the x-direction is not even.'); 
 end
 
  if mod(MESHparams.M, 2) ~= 0; 
     error('The number of intervals in the y-direction is not even.'); 
 end

%% Definition of the SOLVERparmas data structure and its fields

  K =(MESHparams.N-1)*(MESHparams.M-1); % which one?

  SOLVERparams.tol   = 1e-12;      
  SOLVERparams.maxit = K; 
  SOLVERparams.x0    = zeros(K,1);
  SOLVERparams.order = 1; % order of application of the projections in the 
                          % mult. Schwarz iteration operator:
                          % 1: T = Q2Q1, 2: T = Q1Q2.
  
 
end  