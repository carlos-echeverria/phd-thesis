function [ tauu, omegaa ] = TausAndOmegas2( blkA )
%TAUSANDOMEGAS calculates the taus and omegas without refinement.
%
% function call: 
%
%      [ tau, omega ] = TausAndOmegas( blkA )
%
% input:
%
%   blkA: cell array of the M^2 x M^2 matrix A partitioned into blocks of
%         size M x M.
%
% output:
%          
%     tauu: cell array of size M x M where each column corresponds to a
%           refienment step (col 1 - no refinement).
%   omegaa: cell array of size M x M where each column corresponds to a
%           refienment step (col 1 - no refinement)
%
% for the definition of taus and omegas see:
%
%  [ Echeverria, Liesen, Nabben - Block diagonal dominance of matrices
%    revisited: bounds for the norms of inverses and eigenvalue 
%    inclusion sets - Linear Algebra and its Applications - 2018 ] 
%
%   tau_{i,t} =||blkA_{i,i}^{-1}*blkA_{i,i+1}||/(1-||blkA_{i,i}^{-1}*blkA_{i,i-1}||*tau_{i-1,t-1})
% omega_{i,t} =||blkA_{i,i-1}*blkA_{i,i}^{-1}||/(1-||blkA_{i,i+1}*blkA_{i,i}^{-1}||*omega_{i-1,t-1})
%
% Written by Carlos Echeverria. Last edited on May 5, 2017.


M = length(blkA{1,1});

%% build taus and omegas for Theorem 2.6

          tau = cell(M,1);
        omega = cell(M,1);
 
       tau{1} = norm(inv(blkA{1,1})*blkA{1,2});
     omega{1} = 0;
 
 for k = 2:M-1
       tau{k} = norm(inv(blkA{k,k})*blkA{k,k+1})/(1-(norm(inv(blkA{k,k})*blkA{k,k-1})));
     omega{k} = norm(blkA{k,k-1}*inv(blkA{k,k}))/(1-(norm(blkA{k,k+1}*inv(blkA{k,k}))));   
 end
  
     omega{M} = norm(blkA{M,M-1}*inv(blkA{M,M}));
       tau{M} = 0; 

%% build iterative taus and omegas for Theorem 2.7

      tauu = cell(M,M-1);
    omegaa = cell(M,M-1);
 
     for k=1:M-1
          
         tauu{1,k} = norm(blkA{1,1}\blkA{1,2});
         tauu{k,1} = tau{k};
         tauu{M,k} = 0;
         
         omegaa{1,k} = 0; 
         omegaa{k,1} = omega{k};
         omegaa{M,k} = norm(blkA{M,M-1}*inv(blkA{M,M}));
         
     end
     
     for tt=2:M-1
        for ii=2:M-1
            if ii<tt, 
                  tauu{ii,tt} = tauu{ii,tt-1};
            else
                  tauu{ii,tt} = norm(blkA{ii,ii}\blkA{ii,ii+1})/(1-(norm(blkA{ii,ii}\blkA{ii,ii-1})*tauu{ii-1,tt-1}));
            end
          
            if ii>M-tt+1, 
                 omegaa{ii,tt} = omegaa{ii,tt-1};
            else
                  omegaa{ii,tt} = norm(blkA{ii,ii-1}*inv(blkA{ii,ii}))/(1-(norm(blkA{ii,ii+1}*inv(blkA{ii,ii}))*omegaa{ii+1,tt-1}));
            end
        end   
    end
      
end