function [SOLVERparams] = block_Schwarz_2D1L_Kron_get_A( SOLVERparams,MESHparams,EQparams )
%BLOCK_SCHWARZ_2D1L_KRON_GET_LINSYS.
%   Using the equation parameters of a convection-diffusion problem:
%
%  -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma,
%
%   and using the mesh parameters of the discretized domain, this function
%   generates the corresponding coefficient matrix obtained by  using
%   upwind or central difference operators for approximating the first and
%   second order derivatives of the equation. This code uses the
%   Kroenecker product of discrete 1D problems to construct the 2D matrix.
%
%
%   function call:
%
% [SOLVERparams] = block_Schwarz_2D1L_Kron_get_A(SOLVERparams,MESHparams,EQparams)
%
%   input:
%
% SOLVERparams: data structure with solver parameters including rhs 'b'.
%   MESHparams: data structure with enlarged mesh parameters.
%     EQparams: data structure with (at least) the initial equation parameters.
%               (see block_Schwarz_getparams.m for initial params).
%
%   output:
%
%  SOLVERparams: modified data structure. In particular, the coefficient
%                matrix of the limear system is added to the SOLVERparams
%                data structure under the field 'A'.
%                The stored matrix only accounts for the internal points.
%
% Written by Carlos Echeverria on August 8, 2018.
% Last Edited by C. E. on September 25, 2019.


%% Read initial prameters from input data structures

  epsi = EQparams.epsi;
   ph1 = EQparams.ph1;
   ps1 = EQparams.ps1;
   ph2 = EQparams.ph2;
   ps2 = EQparams.ps2;

   Hx  = MESHparams.Hx;
   Hy  = MESHparams.Hy;
   hy  = MESHparams.hy;
%   hx  = MESHparams.hx;


%% Define intermediate sizes for construction of

    N  = MESHparams.N;  % total number of intervals in x-direction
  Nm1  = N-1;           % total number of internal points in x-direction
%    n  = N/2;           % number of local intervals in x-direction
%  nm1  = n-1;           % local number of internal points in x-direction

    M  = MESHparams.M;  % total number of intervals in y-direction
  Mm1  = M-1;           % total number of internal points in y-direction
    m  = M/2;           % number of local intervals in y-direction
  mm1  = m-1;           % local number of internal points in y-direction

  INm1 = eye(Nm1);      % Identity matrix in R^(N-1)
  IMm1 = eye(Mm1);      % Identity matrix in R^(M-1)

% enN  = INm1(:,n);     % n-th cannonical vector in R^(N-1)
  emM  = IMm1(:,m);     % m-th cannonical vector in R^(M-1)


%% Construct Diffusion operator

%  Second derivative operator without mesh scaling:
   T1 = gallery('tridiag',Nm1,1,-2,1);
   T2 = gallery('tridiag',Mm1,1,-2,1);

%  Scaling matrices for second derivative on the regular and 1D Shishkin meshes:
   Dxx = (-1/Hx^2)*diag(ones(Nm1,1));
   Dyy = diag([(-1/Hy^2)*ones(mm1,1);(-1/(Hy*hy));(-1/hy^2)*ones(mm1,1)]);
%  Dyy = (-1/Hy^2)*diag(ones(Mm1,1));
%  Dxx = diag([(-1/Hx^2)*ones(nm1,1);(-1/(Hx*hx));(-1/hx^2)*ones(nm1,1)]);

    Tx = Dxx*T1;
    Ty = Dyy*T2;

%  Rank one correction for middle node in x direction:
   %    gammax = (hx-Hx)/((Hx+hx)*hx*Hx);
   %        vx =  zeros(Nm1,1);
   % vx(n-1,1) = -gammax;
   % vx(n+1,1) =  gammax;

   %        tx = enN*vx';    % need to add if statement for layer in x-dir

%  Rank one correction for middle node in y direction:
      gammay = (hy-Hy)/((Hy+hy)*hy*Hy);
          vy =  zeros(Mm1,1);
   vy(m-1,1) = -gammay;
   vy(m+1,1) =  gammay;

          ty = emM*vy';

%  Second derivative operator on 2D hybrid mesh:
   Diff = epsi*kron(IMm1,Tx) + epsi*kron((Ty+ty),INm1);
%  Diff = epsi*kron(IMm1,(Tx+tx)) + epsi*kron((Ty+ty),INm1);


%% Construct Convection operator

%  First derivative upwind operator on uniform 1D mesh:
   Bup1 = gallery('tridiag',Nm1,-1,1,0);
   Bup2 = gallery('tridiag',Mm1,-1,1,0);
%  Bcd  = gallery('tridiag',M,-1,1,-1);

%  Scaling matrices for first derivative upwind:
   Dx = (1/Hx)*diag(ones(Nm1,1));
   Dy = diag([(1/Hy)*ones(m,1);(1/hy)*ones(mm1,1)]);

   Bupx = Dx*Bup1;
   Bupy = Dy*Bup2;

%  Definition of 'wind' matrices:
   Phi1 = ph1*eye(Nm1);  Psi1 = ps1*eye(Mm1);
   Phi2 = ph2*eye(Nm1);  Psi2 = ps2*eye(Mm1);

%  First derivative operator on 2D hybrid mesh:
   Conv = kron(Psi1,Phi1*(Bupx)) + kron(Psi2*(Bupy),Phi2);


%% Construct Discrete Convection-Diffusion operator

%  Convection-Diffusion operator on 2D Shishkin mesh:
   A_kr = sparse(Conv + Diff);

%% Scaling??
% D = diag([(H/alpha)*ones(n-1,1);...
%           (h*H)/(2*eps);((h*h)/eps)*ones(n-1,1)]); %1D upwind



   SOLVERparams.A = A_kr;

end
