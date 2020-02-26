function [MESHparams] = block_Schwarz_2D1L_mesh( MESHparams,EQparams )
%BLOCK_SCHWARZ_2D1L_MESH generates a piecewise equidistant mesh in the 
% domain Omega = [0,1]x[0,1] for the 2-dimentional problem with maximum of
% 2 boundary layers, one in each coordinate direction.
%
%   Based on the equation parameters of a convection-diffusion problem,
%   this function generates a Shishkin mesh with transition parameters:
%
%               tau_x = min(0.5 ,(2/wx)*epsi*log(N)),
%               tau_y = min(0.5 ,(2/wy)*epsi*log(N)).
%
%   Depending on the components of the wind vector (wx and wy), the 
%   resulting mesh will be refined in the corresponding direction. 
%
%   function call: 
%
%      [MESHparams] = block_Schwarz_2D1L_mesh( MESHparams,EQparams )
%
%   input:
%
%   MESHparams: data structure with the initial mesh parameters:
%          N: even integer with the number of intervals in the x-dir. 
%          M: even integer with the number of intervals in the y-dir.
%
%     EQparams: data structure with the initial equation parameters:

%               see block_Schwarz_getparams.m for initial params.
%
%   output:
%
%   MESHparams: modified data structure. In particular, the following
%               fields are added to the MESHparams data structure:
%               xi: row vector with the x-coordinates of the mesh nodes.
%               yi: row vector with the y-coordinates of the mesh nodes.
%               Hx: coarse mesh parameter in x-direction.
%               hx: fine mesh parameter in x-direction.
%               Hy: coarse mesh parameter in y-direction.
%               hy: fine mesh parameter in y-direction.
%            tau_x: transition parameter in x-direction.
%            tau_y: transition parameter in y-direction.
% 
% Written by Carlos Echeverria on August 8, 2018.
% Edited  by C.E. on September 24, 2019.
% Edited  by C.E. on February  25, 2020.



%% Read off initial prameters from input data structure

  epsi = EQparams.epsi;  
    wx = EQparams.wx;   
    wy = EQparams.wy;
     N = MESHparams.N;
     M = MESHparams.M;
   
     
%% Construct mesh

% x-direction:
    tau_x = min(0.5,2*epsi*log(N)/wx);    %tau_x = 0.5; 
      xi  = zeros(N+1,1);
    
    for i=0:N/2,
        xi(i+1) = 2*i*(1.0-tau_x)/N;
    end

    for i=N/2+1:N,
        xi(i+1) = 1.0 - 2*(N-i)*tau_x/N;
    end
    
    Hx = 2*(1-tau_x)/N;  
    hx = 2*tau_x/N;
    
% y-direction:
    tau_y = min(0.5,2*epsi*log(M)/wy);    %tau_y = 0.5;
      yi  = zeros(M+1,1);
    
    for i=0:M/2,
        yi(i+1) = 2*i*(1.0-tau_y)/M;
    end
    
    for i=M/2+1:M,
        yi(i+1) = 1.0 - 2*(M-i)*tau_y/M;
    end
    
    Hy = 2*(1-tau_y)/M;  
    hy = 2*tau_y/M;
       
    
%% Save prameters of constructed mesh to update initial data structure

    MESHparams.xi = xi;
    MESHparams.yi = yi;
    MESHparams.Hx = Hx;
    MESHparams.hx = hx;
    MESHparams.Hy = Hy;
    MESHparams.hy = hy;
    MESHparams.tau_x = tau_x;
    MESHparams.tau_y = tau_y;
    
end