function [ SOLUTIONparams ] = block_Schwarz_error_bounds_model_problems_inexact(EQparams, MESHparams, SOLUTIONparams )
%BLOCK_SCHWARZ_ERROR_BOUNDS computes the error generated by the
% multiplicative Schwarz method at each iteration as well as the error
% bounds to the Schwarz iterations.
%
%   function call:
%
% [e_boundT12, a_boundT21] = ...
%          block_Schwarz_error_bounds2(EQparams,MESHparams,SOLUTIONparams)
%
%   input:
%
%        EQparams: data structure array with equation parameters.
%      MESHparams: data structure with mesh parameters.
%  SOLUTIONparams: data structure with solution parameters.
%
%   output:
%
%      err_boundT12: vector with the values of the theoretical error bound
%                    for row and column block diagonal dominat matrices
%                    for the mult. Schwarz method with ordering T = Q2Q1.
%      err_boundT21: vector with the values of the theoretical error bound
%                    for row and column block diagonal dominat matrices
%                    for the mult. Schwarz method with ordering T = Q1Q2.
%
% Written by Carlos Echeverria on January 19, 2019
% Last Edited by C. E. on September 24, 2019.



%%  Compute norms(Tij) and convergence factor for Theorem 4.4

if EQparams.problem == 1

  %  Define needed entries of the matrix for convection diffusion bound:
     eH = -EQparams.epsi/MESHparams.Hy^2;
     dH = (-EQparams.epsi/MESHparams.Hy^2)-(EQparams.wy/MESHparams.Hy);

     %          rho = abs(eH/dH);
     % norm_T12_inf = abs(eH/dH);

                rho = abs(EQparams.epsi/(EQparams.epsi-MESHparams.Hy));
       norm_T12_inf = abs(EQparams.epsi/(EQparams.epsi-MESHparams.Hy));

       norm_T21_inf = 1;

elseif EQparams.problem == 2

  %  Define needed quantities for Laplace bound:
     l1 = 4-2*cos(pi/(MESHparams.N)); % eigenvalue of T=tridiag(-1,2,-1)

          rho = (1/(l1^2-l1-1))^2;

       norm_T12_inf = 1/l1;
       norm_T21_inf = 1/l1;

end

%%  Compute bounds of Theorem 4.4

% bounds for T12 (vectorized):
iter=1:length(SOLUTIONparams.err_s1)-2;
err_boundT12 =[1, norm_T12_inf, (rho.^iter)*norm_T12_inf ];

SOLUTIONparams.err_bound_s2 = err_boundT12;
SOLUTIONparams.rho = rho;

end
