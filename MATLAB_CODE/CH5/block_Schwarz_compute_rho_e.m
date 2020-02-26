function [A,B,C,AH,BH,CH,Ah,Bh,Ch,rho_e_inf] = block_Schwarz_compute_rho_e( A, MESHparams )
%BLOCK_SCHWARZ_COMPUTE_RHO_E computes the "exact" convergence factor of the
% multiplicative Schwarz method given in eq (5.12) of CH.5 of the 
% PH.D. thesis:
%
%  [ Echeverria - Iterative solution of discretized convection-diffusion 
%                 problems - Technische Universitaet Berlin - 2020 ]
%
% The factor is calculated by computing the norms of the relevant blocks
%  of the iteration matrix of the method.
%
%   function call:
%
% [A,B,C,AH,BH,CH,Ah,Bh,Ch,rho_e_inf] = block_Schwarz_compute_rho_e(A,MESHparams)
%
%   input:
%
%              A: block matrix obtained from discretizaton.
%     MESHparams: data structure with mesh parameters.
%
%   output:
%
%       rho_e_inf: convergence factor calculated using the infinity norm
%
% Written by Carlos Echeverria on November 5, 2019.


Nm1 = MESHparams.N-1;
Mm1 = MESHparams.M-1;
 M2 = MESHparams.M/2;
 N2 = MESHparams.N/2;

%% Construct iteration matrix of the multiplicative Schwarz method

% Notation conflict: In the physical problem we have that each local
%   subdomain consists of M/2 blocks of size (N-1)x(N-1). To be consistent
%   with the notatoion of the restriction operators in Ch5 of the thesis
%   we now create local variables where (N-1) is referred to as N and (M/2)
%   is referred to as (m+1):

 N = Nm1;     n = M2;       m = n-1;    M = Mm1;

% Extract the submatrices AH and Ah:
   AH = A(1:N*m,1:N*m);
   Ah = A(N*n+1:M*N,N*n+1:M*N);

%%  Invert block matrix 'A' and also the submatrices 'AH' and 'Ah':

      Z = A\eye(N*M);
     ZH = AH\eye(N*m);
     Zh = Ah\eye(N*m);

  Zdev  = norm(Z*A-eye(N*M));     % check if Z  fullfils   A*Z=I
 ZHdev  = norm(ZH*AH-eye(N*m));   % check if ZH fullfils AH*ZH=I
 Zhdev  = norm(Zh*Ah-eye(N*m));   % check if Zh fullfils Ah*Zh=I


%% partition matrices 'A' and 'Z' into blocks and store them in cells:

   [  blkA ] = partition_block_matrix(  A, N, M);
   [  blkZ ] = partition_block_matrix(  Z, N, M);
   [ blkAH ] = partition_block_matrix( AH, N, m);
   [ blkZH ] = partition_block_matrix( ZH, N, m);
   [ blkAh ] = partition_block_matrix( Ah, N, m);
   [ blkZh ] = partition_block_matrix( Zh, N, m);

%% Extract blocks of the matrices

     AH = blkAH{1,1};    BH = blkAH{1,2};     CH = blkAH{2,1};
     Ah = blkAh{1,1};    Bh = blkAh{1,2};     Ch = blkAh{2,1};
     A  = blkA{n,n};     B  = blkA{n,n+1};    C  = blkA{n,n-1};
   ZHmm = blkZH{m,m};  Zh11 = blkZh{1,1};

% % note: this snippet of code only works in case of block tridiagonal Toeplitz matrices

Pi1=(A-C*ZHmm*BH)\B;
Pi2=(A-B*Zh11*Ch)\C;

rho_e_inf=norm(Zh11*Ch*Pi2*ZHmm*BH*Pi1, inf);
rho_e_two=norm(Zh11*Ch*Pi2*ZHmm*BH*Pi1, 2);

rho_ex =norm(ZHmm*BH*Pi1*Zh11*Ch*Pi2, inf);
rho_ex_two=norm(ZHmm*BH*Pi1*Zh11*Ch*Pi2, 2);

% fprintf('rho_12_inf:%5.15e, rho_12_two:%5.15e, rho_21_inf:%5.15e, rho_21_two:%5.15e \n', rho_e_inf, rho_e_two, rho_ex, rho_ex_two)


end
