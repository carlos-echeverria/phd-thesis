function [ FVVa, FVVb, ELNN1a, ELNN1b, ELNN2a, ELNN2b ] = BlockDiagDomEval( blkA, I, x, y )
%BDIDOEVAL Summary of this function goes here
%
%
%   function call:
%
%  [ FVV, ELNN1, ELNN2 ] = BDiDoEval( blkA, I, x, y )
%
%      input:
%
%    blkA: cell array of the M^2 x M^2 matrix A partitioned into blocks of
%           size M x M.
%       I: Identity of size M x M
%       x: one directional coordinates of x-direction
%       y: one directional coordinates of y-direction
%
%      output:
%
%     FVV: the function of the complex plane evaluated at the points z=x+iy
%     LN1: the function of the complex plane evaluated at the points z=x+iy
%     LN2: the function of the complex plane evaluated at the points z=x+iy
%
% Written by Carlos Echeverria. Last edited on Nov. 25, 2017.


%% define functions  FV and ELN1 & ELN2

    FVa = @(x,y) norm(blkA{1,2})*norm((blkA{1,1}-(x+1i*y).*I)\I);
    FVb = @(x,y) norm(blkA{2,1})*norm((blkA{2,2}-(x+1i*y).*I)\I);

  ELN1a = @(x,y) norm((blkA{1,1}-(x+1i*y).*I)\blkA{1,2});
  ELN1b = @(x,y) norm((blkA{2,2}-(x+1i*y).*I)\blkA{2,1});
  
  ELN2a = @(x,y) norm(blkA{1,2}/(blkA{1,1}-(x+1i*y).*I));
  ELN2b = @(x,y) norm(blkA{2,1}/(blkA{2,2}-(x+1i*y).*I));


   FVVa = zeros(size(x,2),size(y,2));
   FVVb = zeros(size(x,2),size(y,2));
 ELNN1a = zeros(size(x,2),size(y,2));
 ELNN1b = zeros(size(x,2),size(y,2));
 ELNN2a = zeros(size(x,2),size(y,2));
 ELNN2b = zeros(size(x,2),size(y,2));

for i=1:size(x,2)
    for j=1:size(y,2)

        FVVa(j,i)= FVa(x(i),y(j));
        FVVb(j,i)= FVb(x(i),y(j));

      ELNN1a(j,i)= ELN1a(x(i),y(j));
      ELNN1b(j,i)= ELN1b(x(i),y(j));

      ELNN2a(j,i)= ELN2a(x(i),y(j));
      ELNN2b(j,i)= ELN2b(x(i),y(j));

    end
end

end
