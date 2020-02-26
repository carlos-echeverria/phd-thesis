function [ blkA ] = PartitionMatrix( A )
%PARTITIONMATRIX uses the cell data type to partition an M^2 x M^2 matrix
% 'A' into cell blocks of size M x M. Block ij can be accesed by blkA{i,j}.
%
% function call: 
%
%   [ blkA ] = BlockTridiagBounds( A )
%
% input:
%
%        A: matrix of size M^2 x M^2.
%
% output:
%          
%     blkA: M x M array of cells, each one of size M x M.  
%
% Written by Carlos Echeverria. Last edited on Nov. 25, 2017.

M = sqrt(size(A,1));

blkA = cell(M,M);
for i=1:M
    for j=1:M
        blkA{i,j}=full(A(((i-1)*M+1):(i*M),((j-1)*M+1):j*M));
        %blockA(i,j)=cell2mat(blkA())
    end
end

end

