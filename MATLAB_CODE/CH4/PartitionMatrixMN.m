function [ blkA ] = PartitionMatrixMN( A, M, N )
%PARTITIONMATRIXMN uses the cell data type to partition an MN x MN matrix
% 'A' into N cell blocks of size M x M. Block ij can be accesed by
%  blkA{i,j}.
%
% function call:
%
%   [ blkA ] = BlockTridiagBoundsMN( A, M, N )
%
% input:
%
%   A: blok matrix of size M*N x M*N.
%
% output:
%
%   blkA: N x N array of cells, each one of size M x M.
%
% Written by Carlos Echeverria on October 24, 2019.


blkA = cell(N,N);

for i=1:N
    for j=1:N
        blkA{i,j}=full(A(((i-1)*M+1):(i*M),((j-1)*M+1):j*M));
        %blockA(i,j)=cell2mat(blkA())
    end
end

end
