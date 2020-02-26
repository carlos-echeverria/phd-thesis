function [ Z, U, V, X, Y ] = inverseUVXY( blkA )
%INVERSEUVXY calculates the inverse of a block tridiagonal matrix blkA 
% by constructing sequences of matrices U, V, X, and Y, following:
%
%[Ikebe Y. - On Inverses of Hessenberg Matrices, Linear Algebra Appl. 1971]
%
% function call: 
%
%             [ Z ] = inverseUVXY( blkA )
%
% input:
%
%   blkA: cell array of the M^2 x M^2 matrix A partitioned into blocks of
%         size M x M. (see PartitionMatrix.m)
%
% output:
%          
%     Z: M x M cell array containing the inverse of A partitioned into
%           blocks of size M x M.  
%
% Written by Carlos Echeverria. Last edited on April 30, 2018.

M = length(blkA{1,1});

%% generate sequence of matrices 'U'

U{1} = eye(M);

U{2} = -blkA{1,2}\blkA{1,1}*U{1}; 

for i=3:M
    MAT3 = zeros(M,M);
        for k=1:i-1
            MAT3 = MAT3 + blkA{i-1,k}*U{k};
        end
      U{i} = -blkA{i-1,i}\MAT3; 
end
      

%% generate sequence of matrices 'V'

SUM2 = zeros(M,M);
 for k=1:M
     SUM2 = SUM2+blkA{M,k}*U{k};
 end
V{M} = inv(SUM2); 

V{M-1} = -V{M}*blkA{M,M}/blkA{M-1,M};

for i=M-2:-1:1                              
    MAT4 = zeros(M,M);                       
        for k=i+1:M
            MAT4 = MAT4+V{k}*blkA{k,i+1};        
        end
        %pause
      V{i} = -MAT4/blkA{i,i+1};
end
 

%% generate sequence of matrices 'X'

X{1} = eye(M);

X{2} = -X{1}*blkA{1,1}/blkA{2,1};

for i=3:M
    MAT = zeros(M,M);
        for k=1:i-1
            MAT = MAT+X{k}*blkA{k,i-1};
        end
      X{i} = -MAT/blkA{i,i-1};
end


%% generate sequence of matrices 'Y'  

SUM = zeros(M,M);
 for k=1:M
     SUM = SUM + X{k}*blkA{k,M};
 end
Y{M} = inv(SUM);

Y{M-1} = -blkA{M,M-1}\blkA{M,M}*Y{M};

for i=M-2:-1:1
    MAT2 = zeros(M,M);
        for k=i+1:M
            MAT2 = MAT2 + blkA{i+1,k}*Y{k};
        end
      Y{i} = -blkA{i+1,i}\MAT2;
end


%% compute matrix 'Z' which is the inverse of 'A'

Z = cell(M,M);
for i=1:M
    for j=1:M
        if i<=j, Z{i,j} = U{i}*V{j};end
        if i>j,  Z{i,j} = Y{i}*X{j};end
    end
end

end