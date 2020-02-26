function [BDD, BDDFV] = BlockDiagDomCheck( blkA )
%BLOCKDIAGDOMCHECK checks if the cell array blkA is block diagonally 
% dominant with respect to both FV's and ELN's definitions.
%
% function call: 
%
%           [ blkA ] = BlockDiagDomCheck( blkA )
%
% input:
%
%   blkA: cell array of the M^2 x M^2 matrix A partitioned into blocks of
%         size M x M.
%
% output:
%          
%    BDD: vector with the quantity \sum_{i\neq j}||A_{ii}^{-1}A_{ij}||
%                for each block row i=1 to M.
%
%  BDDFV: vector with the quantity \sum_{i\neq j}||A_{ij}||?||A_{ii}^{-1}||
%                for each block row  i=1 to M.
%
% function returns an error if the matrix is not block diagonally domiant.
%
% Written by Carlos Echeverria. Last edited on April. 30, 2017.

M = length(blkA{1,1});

  BDD = zeros(1,M);
BDDFV = zeros(1,M);

  BDD(1) = norm(blkA{1,1}\blkA{1,2});
BDDFV(1) = norm(inv(blkA{1,1}))*norm(blkA{1,2});

if BDD(1)>1, error('Matrix A is not row block diagonally dominant.'); end

for k=2:M-1
      BDD(k) = norm(blkA{k,k}\blkA{k,k-1}) + norm(blkA{k,k}\blkA{k,k+1});
    BDDFV(k) = norm(blkA{k,k-1})*norm(inv(blkA{k,k})) + norm(blkA{k,k+1})*norm(inv(blkA{k,k}));
    if BDD(k)>1, error('Matrix A is not row block diagonally dominant.'); end
end

  BDD(M) = norm(blkA{M,M}\blkA{M,M-1});
BDDFV(M) = norm(inv(blkA{M,M}))*norm(blkA{M,M-1});

if BDD(M)>1, error('Matrix A is not row block diagonally dominant.'); end

end

