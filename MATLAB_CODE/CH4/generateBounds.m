function [ UBound, LBound, upbound, lowbound, entries, diagentries ] = generateBounds(blkZ, blkA, tauu, omegaa, t )
%GENERATEBOUNDS. Using the computed taus and omegas, this function generates
% the bounds on the norms of the off-diagonal as well as the diagonal 
% block elements of the block tridiagonal matrix A for the specific 
% iterative refinement step t.
%
%   function call: 
%
% [ UBound, LBound, upbound, lowbound, entries, diagentries ] = generateBounds(blkZ, blkA, tauu, omegaa, t )
%
%    input:
%
%    blkZ: cell array of the M^2 x M^2 matrix Z partitioned into blocks of
%          size M x M (Z is the inverse of A).
%    blkA: cell array of the M^2 x M^2 matrix A partitioned into blocks of
%          size M x M.
%    tauu: M x M cell array containing the computed taus for all refienement
%          steps.
%  omegaa: M x M cell array containing the computed omegas for all 
%          refienement steps.
%       t: number between 1 and M-1 specifying the level of refinement.
%
%   output:
%          
%      UBound: M x M matrix with the upper bounds on the norms of the blocks
%              of Z in its entries for the specified refinement step.
%      LBound: M x M diagonal matrix with the lower bounds on the norms of 
%              the blocks of Z in its entries for the specified refinement 
%              step.
%    lowbound: lower triangular matrix with the upper bounds on the lower 
%              triangular off-diagonal blocks of Z.
%     upbound: upper triangular matrix with the upper bounds on the upper 
%              triangular off-diagonal blocks of Z.
%     entries: M x M matrix with the exact norms of the blocks of Z in its
%              entries.
% diagentries: 1 x M vector with the exact norms of the diagonal blocks
%              of Z.
%
% Written by Carlos Echeverria. Last edited on Nov. 25, 2017.


M = length(blkA{1,1});


%% calculate upper bounds for the off-diagonal entries:

 entries = zeros(M,M);
  ubound = zeros(M,M);
  lbound = zeros(M,M);
  
 for i=1:M
     for j=1:M
         
        entries(i,j) = norm(blkZ{i,j});
        
        
        if i<j 
             efftau = 1;
             for k = i:j-1
                  efftau = efftau*tauu{k,t};
             end
            ubound(i,j) = norm(blkZ{j,j})*efftau;
        end
 
        
        if i>j 
           effomega = 1;
           for k = j+1:i
               effomega = effomega*omegaa{k,t};
           end
            lbound(i,j) = norm(blkZ{j,j})*effomega;
        end
        
        
     end
 end

%% calculate upper and lower bounds for the diagonal entries:

      lowbound = cell(M);
       upbound = cell(M);
   diagentries = zeros(M,1);
 
 diagentries(1) = norm(blkZ{1,1});
    lowbound{1} = 1/(norm(blkA{1,1})+(omegaa{2,t}*norm(blkA{1,2})));
     upbound{1} = 1/((1/(norm(inv(blkA{1,1}))))-(omegaa{2,t}*norm(blkA{1,2})));
 
 for i=2:M-1
     
     diagentries(i) = norm(blkZ{i,i});
     
               a = norm(blkA{i,i});
              aa = 1/norm(inv(blkA{i,i}));
               b = tauu{i-1,t}*norm(blkA{i,i-1});
               c = omegaa{i+1,t}*norm(blkA{i,i+1});
     lowbound{i} = 1/(a+b+c);
      upbound{i} = 1/(aa-b-c);
     
 end
 
 diagentries(M) = norm(blkZ{M,M});
    lowbound{M} = 1/(norm(blkA{M,M})+(tauu{M-1,t}*norm(blkA{M,M-1})));
     upbound{M} = 1/(1/(norm(inv(blkA{M,M})))-(tauu{M-1,t}*norm(blkA{M,M-1})));

     
lowbound = cell2mat(lowbound);
 upbound = cell2mat(upbound);
 
  UBound = ubound + lbound + diag(upbound);
  LBound = diag(lowbound);
  
end

