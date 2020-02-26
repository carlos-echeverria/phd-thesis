function [ RBDD_boundT12, RBDD_boundT21, RCBDD_boundT12, RCBDD_boundT21 ] = block_Schwarz_error_bounds(blkA, blkAH, blkAh, err_s1, err_s2 )
%BLOCK_SCHWARZ_ERROR_BOUNDS computes the theoretical error bounds for the  
% Multiplicative Schwarz method given in: 
%
% [ Echeverria, Liesen, Tichy - Analysis of the multiplicative Schwarz 
%   method for matrices with a special block structure - 2019 
%   manuscript version: September 6, 2019 ]
%
% The bounds are calculated by first computing intermediate quantities like
% convergence factors and error norms of the iteration matrix of the method.
%
%   function call:
%
%    [ e_bound1, a_bound1, e_bound2, a_bound2 ] = ...
%          block_Schwarz_error_bounds(blkA, blkAH, blkAh, err_s1, err_s2 )
%
%   input: 
%
%           blkA: cell array of the matrix A  partitioned into blocks.
%          blkAh: cell array of the matrix Ah partitioned into blocks.
%          blkAH: cell array of the matrix AH partitioned into blocks.
%         err_s1: infinity norm of errors at each iteration for T = Q2Q1.
%         err_s2: infinity norm of errors at each iteration for T = Q1Q2.
%
%   output:
%
%       RBDD_boundT12: vector with the bound for row-block diagonal dominat matrices T = Q2Q1.                
%       RBDD_boundT21: vector with the bound for row-block diagonal dominat matrices T = Q1Q2.
%      RCBDD_boundT12: vector with the bound for row and column-block diagonal dominat matrices T=Q2Q1.                
%      RCBDD_boundT21: vector with the bound for row and column-block diagonal dominat matrices T=Q1Q2.  
%
%
% Written by Carlos Echeverria on August 7, 2018.
% Last Edited by C. E. on September 24, 2019.

 m = length(blkAh);
 n = m + 1;
 
%% Extract blocks of the matrices

     AH = blkAH{1,1};    BH = blkAH{1,2};     CH = blkAH{2,1};
     Ah = blkAh{1,1};    Bh = blkAh{1,2};     Ch = blkAh{2,1};
     A  = blkA{n,n};     B  = blkA{n,n+1};    C  = blkA{n,n-1};
     
% note: this snippet of code only works in case of block tridiagonal Toeplitz matrices
     
%% Calculate needed quantities for bounds of Theorem 4.3 (manuscript version: 06.09.19)

      eta_h = (norm(Ah\Ch))/(1-norm(Ah\Bh)); % Eq.(4.5)
      eta_H = (norm(AH\BH))/(1-norm(AH\CH)); % Eq.(4.6)
      
  eta_h_inf = (norm(Ah\Ch,inf))/(1-norm(Ah\Bh,inf)); % Eq.(4.11)
  eta_H_inf = (norm(AH\BH,inf))/(1-norm(AH\CH,inf)); % Eq.(4.11)
      
        rho = ((eta_h*norm(A\C))/(1-(eta_h*norm(A\B))))*((eta_H*norm(A\B))/(1-(eta_H*norm(A\C)))); % Thm. 4.3
   
   norm_T12_inf = (eta_H_inf*norm(A\B,inf))/(1-(eta_H_inf*norm(A\C,inf))); % Thm. 4.3
   norm_T21_inf = (eta_h_inf*norm(A\C,inf))/(1-(eta_h_inf*norm(A\B,inf))); % Thm. 4.3
       

%% Compute bounds of Theorem 4.3 for block-row diagonal dominat matrices 

% bounds for T12 (vectorized):
k = 1:length(err_s1)-2; % correct length of iteration steps
RBDD_boundT12 = [1, norm_T12_inf, (rho.^k)*norm_T12_inf]; 


% bounds for T21 (vectorized):
k = 1:length(err_s2)-2;
RBDD_boundT21 = [1, norm_T21_inf, (rho.^k)*norm_T21_inf];


%%  Calculate needed quantities for bounds of Theorem 4.4 (manuscript version: 06.09.19)
      
   eta_h_inf_2 = (norm(inv(Ah), inf)*norm(Ch, inf))/(1-norm(Ch/Ah, inf)); % Eq. (4.13)
   eta_H_inf_2 = (norm(inv(AH), inf)*norm(BH, inf))/(1-norm(BH/AH, inf)); % Eq. (4.14)
    
 eta_h_min_inf = min(eta_h_inf, eta_h_inf_2); % Eq. (4.13)
 eta_H_min_inf = min(eta_H_inf, eta_H_inf_2) ;% Eq. (4.14)
    
       eta_h_2 = (norm(inv(Ah))*norm(Ch))/(1-norm(Ch/Ah)); % Eq. (4.13)
       eta_H_2 = (norm(inv(AH))*norm(BH))/(1-norm(BH/AH)); % Eq. (4.14)
    
     eta_h_min = min(eta_h, eta_h_2); % Eq. (4.13)
     eta_H_min = min(eta_H, eta_H_2); % Eq. (4.14)
    
         rho_2 = ((eta_h_min*norm(A\C))/(1-(eta_h_min*norm(A\B)))) * ((eta_H_min*norm(A\B))/(1-(eta_H_min*norm(A\C)))); % Thm. 4.3
    
   norm_T12_min_inf = (eta_H_min_inf*norm(A\B,inf))/(1-(eta_H_min_inf*norm(A\C,inf))); % Thm. 4.3
   norm_T21_min_inf = (eta_h_min_inf*norm(A\C,inf))/(1-(eta_h_min_inf*norm(A\B,inf))); % Thm. 4.3


 %%  Calculate bounds of Theorem 4.4  for block-row and block-column diagonal dominant matrices
 
% bounds for T12 (vectorized):
k = 1:length(err_s1)-2;
RCBDD_boundT12 = [1, norm_T12_min_inf, (rho_2.^k)*norm_T12_min_inf]; 

% bounds for T21 (vectorized):
k = 1:length(err_s2)-2;
RCBDD_boundT21 = [1, norm_T21_min_inf, (rho_2.^k)*norm_T21_min_inf]; 


end