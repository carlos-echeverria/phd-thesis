function [commutativityAB, commutativityAC] = CommutativityCheck( blkA )
%COMMUTATIVITYCHECK checks if the blocks of the cell array blkA 
% commute with each other.
%
% function call: 
%
%           [ blkA ] = CommutativityCheck( blkA )
%
% input:
%
%   blkA: cell array of the M^2 x M^2 matrix A partitioned into blocks of
%         size M x M.
%
% output:
%          
%  commutativityAB: flag vector of size 1 x M which checks the 
%                   commutativity of the blocks A_i and B_i for each i.
%                   the entry of the vector is 1 if the blocks commute and 
%                   0 if the blocks do not comute.
%
%  commutativityAC: flag vector of size 1 x M which checks the 
%                   commutativity of the blocks A_i and C_{i-1} for each i.
%                   the entry of the vector is 1 if the blocks commute and 
%                   0 if the blocks do not comute.
%
% Written by Carlos Echeverria. Last edited on May 6, 2018.

M = length(blkA{1,1});

commutativityAB = zeros(1,M);
commutativityAC = zeros(1,M);

commutativityAB(1) = isequal(blkA{1,2}*blkA{1,1},blkA{1,1}*blkA{1,2});
     commutativityAC(1) = 1;


for i=2:M-1
     commutativityAB(i) = isequal(blkA{i,i+1}*blkA{i,i},blkA{i,i}*blkA{i,i+1});
     commutativityAC(i) = isequal(blkA{i,i-1}*blkA{i,i},blkA{i,i}*blkA{i,i-1});

end

    commutativityAB(M) = 1;
    commutativityAC(M) = isequal(blkA{M,M-1}*blkA{M,M},blkA{M,M}*blkA{M,M-1});

end