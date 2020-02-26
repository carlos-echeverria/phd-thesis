function [ SOLVERparams ] = block_Schwarz_2D1L_get_b(SOLVERparams,MESHparams, EQparams )
%BLOCK_SCHWARZ_2D1L_GET_B. 
%   Based on the equation parameters of a 2D1L convection-diffusion problem
%   this function evaluates the right hand side of the PDE, 'f', on a given 
%   set of points (xi,yi). The boundary condition function, 'g', is then
%   used to incorporate the contribution of the boundary conditions to a 
%   desired set of points.
%
%   function call: 
%
% [SOLVERparams] = block_Schwarz_2D1L_get_b(SOLVERparams,MESHparams,EQparams)
%
%   input:
%
% SOLVERparams: data structure with (at least) the initial solver parameters.
%   MESHparams: data structure with updated mesh parameters.
%     EQparams: data structure with (at least) the initial equation parameters.
%               (see block_Schwarz_getparams.m for initial params).
%
%   output:
%
%  SOLVERparams: modified data structure. In particular, the right hand
%                side vector of the limear system is added to the 
%                SOLVERparams data structure under the field 'b'. 
%                The stored vector does not include the boundary points.
%
% Written by Carlos Echeverria on August 8, 2018.
% Edited  by C.E. on September 25, 2019.
% Edited  by C.E. on February  25, 2020.



%% Convection-Diffusion Problem (Nonzero boundary function)


if EQparams.problem == 1
    
 % Read equation parameters
    
      epsi = EQparams.epsi; 
        wx = EQparams.wx; 
        wy = EQparams.wy;
         f = EQparams.f;
         g = EQparams.g;
         
 % Read mesh parameters
    
        N  = MESHparams.N;
        M  = MESHparams.M;
        xi = MESHparams.xi;
        yi = MESHparams.yi;
        Hx = MESHparams.Hx;
        Hy = MESHparams.Hy;
        hy = MESHparams.hy;
    
 % Create grid from meshpoints
   
     [X,Y] = meshgrid(xi,yi); 

 %  Apply contribution of the nonzero boundary conditions to discrete RHS
      
      b = f(X',Y'); 

      %  contribution of south boundary 
      b(2:N,2) = b(2:N,2) + g(xi(2:N),0 )*epsi/Hy^2 + g(xi(2:N),0)*wy/Hy;
      
      %  contribution of north boundary 
      b(2:N,M) = b(2:N,M) + g(xi(2:N),1 )*epsi/hy^2;       

      %  contribution of west boundary
      b(2,2:M) = b(2,2:M) + g(0,yi(2:M)')*epsi/Hx^2 + g(0,yi(2:M)')*wx/Hx; 
      
      %  contribution of east boundary
      b(N,2:M) = b(N,2:M) + g(1,yi(2:M)')*epsi/Hx^2;      
      
      
      b = b(2:N,2:M);                    % remove boundary points
      b = reshape(b,[(N-1)*(M-1),1]);    % reshape to column vector
      SOLVERparams.b = b;                % store value
end



%% Laplace Problem 

if EQparams.problem == 2
        
       N = MESHparams.N;  
       M = MESHparams.M;    
      Hx = MESHparams.Hx;
      Hy = MESHparams.Hy;
      hy = MESHparams.hy;

    epsi = EQparams.epsi;      
       f = EQparams.f;
       g = EQparams.g;

      xi = MESHparams.xi;
      yi = MESHparams.yi;
    
   [X,Y] = meshgrid(xi,yi); 

      
%  Apply contribution of the boundary conditions to the discrete RHS
           
      b = f(X',Y');
            
      b(2:N,2) = b(2:N,2) + g(xi(2:N), 0)*epsi/Hy^2; % south boundary 
      b(2:N,M) = b(2:N,M) + g(xi(2:N), 1)*epsi/hy^2; % north boundary      
      b(2,2:M) = b(2,2:M) + g(0,yi(2:M)')*epsi/Hx^2; % west  boundary
      b(N,2:M) = b(N,2:M) + g(1,yi(2:M)')*epsi/Hx^2; % east  boundary
      
      b = b(2:N,2:M);                    % remove boundary points
      b = reshape(b,[(N-1)*(M-1),1]);    % reshape to column vector 
      SOLVERparams.b = b;                % store value
     
end


end