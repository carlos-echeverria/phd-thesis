function [AH,Ah,Z,ZH,Zh,blkA,blkAh,blkAH,blkZ,blkZH,blkZh,Zdev,ZHdev,Zhdev]  = block_Schwarz_extract_and_partition(A, MESHparams)
%BLOCK_SCHWARZ_EXTRACT_AND_PARTITION. Given the block matrix: 
%
%                       [  AH   | [0 BH]' |        ]
%                       [ [0 C] |    A    | [B 0]  ]
%                       [       | [CH 0]' |   Ah   ]
%
% This function extracts the matrices AH, Ah and partitions them into
% blocks with size relative to the parameters of the mesh (size of domain). 
% The inverses of the matrices A, Ah and AH are computed and also
% partitioned into blocks. The matrices are then returned both as cell
% arrays and as MATLAB matrices.
%
%
%   function call: 
%
% [AH,Ah,Z,ZH,Zh,blkA,blkAh,blkAH,blkZ,blkZH,blkZh,Zdev]  = ...
%                       block_Schwarz_extract_and_partition(A, MESHparams)
%
%     input:  
%
%              A: Matrix to be partitioned into blocks.
%     MESHparams: data structure with mesh parameters.
%
%    output:
%
%      AH: matrix of the inverse of AH
%      Ah: matrix of the inverse of Ah
%       Z: matrix of the inverse of A
%      ZH: cell array of the inverse of AH
%      Zh: cell array of the inverse of Ah
%    blkA: cell array of the matrix A partitioned into blocks.
%   blkAh: cell array of the matrix Ah partitioned into blocks.
%   blkAH: cell array of the matrix AH partitioned into blocks.
%    blkZ: cell array of the inverse of A partitioned into blocks.
%   blkZH: cell array of the inverse of AH partitioned into blocks.
%   blkZh: cell array of the inverse of Ah partitioned into blocks.
%    Zdev: measure of the exactness of the calculated inverse of A.
% 
% Written by Carlos Echeverria on September 13, 2018.
% Last edited by C.E. on September 18, 2019.


  M = MESHparams.M; 
  m = MESHparams.m; 
  n = MESHparams.n; 

%% Extract the submatrices AH and Ah:
 
   AH = A(1:M*m,1:M*m);
   Ah = A(M*n+1:M*M,M*n+1:M*M);
 
   
%%  Invert block matrix 'A' and also the submatrices 'AH' and 'Ah':
            
      Z = A\eye(M*M);          
     ZH = AH\eye(M*m);
     Zh = Ah\eye(M*m);
    
  Zdev  = norm(Z*A-eye(M*M));     % check if Z  fullfils   A*Z=I
 ZHdev  = norm(ZH*AH-eye(M*m));   % check if ZH fullfils AH*ZH=I
 Zhdev  = norm(Zh*Ah-eye(M*m));   % check if Zh fullfils Ah*Zh=I


 
%% partition matrices 'A' and 'Z' into blocks and store them in cells:

   [  blkA ] = partition_block_matrix(  A, M, M);
   [  blkZ ] = partition_block_matrix(  Z, M, M);
   [ blkAH ] = partition_block_matrix( AH, M, m);
   [ blkZH ] = partition_block_matrix( ZH, M, m);
   [ blkAh ] = partition_block_matrix( Ah, M, m);
   [ blkZh ] = partition_block_matrix( Zh, M, m);
   
   
end