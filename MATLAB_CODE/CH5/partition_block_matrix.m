function [ blkA ] = partition_block_matrix( A, M, N )
%PARTITION_BLOCK_MATRIX uses the cell data type to partition an NM x NM 
%  matrix 'A' into cell blocks of size M x M. 
%  The ij-th Block can be accesed by the command: blkA{i,j}.
%
% function call: 
%
%      [ blkA ] = partition_block_matrix( A, M, N )
%
% input:
%
%        A: block matrix of size MN x MN.
%
% output:
%          
%     blkA: N x N cell array, each cell of size M x M containing the blocks
%           of A.
%
% Written by Carlos Echeverria on Nov. 25, 2017.
% Edited by C. E. on Sept. 12, 2018.


blkA = cell(N,N);

for i=1:N
    for j=1:N
        blkA{i,j}=full(A(((i-1)*M+1):(i*M),((j-1)*M+1):j*M));
        %blockA(i,j)=cell2mat(blkA())
    end
end


end