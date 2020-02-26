function [SOLUTIONparams] = block_Multiplicative_Schwarz(SOLVERparams, MESHparams, SOLUTIONparams)
% BLOCK_MULTIPLICATIVE_SCHWARZ is an implementation of the 
% mutiplicative Schwarz method for solving two subdomains.
%
%   function call: 
%
% [SOLUTIONparams] = block_Multiplicative_Schwarz(MESHparams,SOLVERparams,SOLUTIONparams)
%
%   input:
%
%     MESHparams: data structure with mesh parameters.
%   SOLVERparams: data structure with solver parameters.
% SOLUTIONparams: data structure with solution parameters.
%                 (see block_Schwarz_getparams.m for initial params).
%
%
%   output:
%
% SOLUTIONparams: data structure with the following solution parameters.
%                 x_s: computed approximate solution of linear system. 
%               res_s: vector with the values of the residulas at each iteration.
%               err_s: infinity norms of the errors at each iteration.
%                   c: correction vector.
%                   T: Multiplicative Scwarz iteration matrix.
%                  Pi: Projection  operator to the i-th domain.
%                  Ri: Restriction operator to the i-th domain.
%
% Written by Carlos Echeverria on January 31, 2016.
% Edited  by Joerg Liesen on September 12, 2016.
% Edited  by C.E. on September 24, 2019.
% Edited  by C.E. on February  25, 2020.



%% Extract needed values from data structures

  tol = SOLVERparams.tol;      
maxit = SOLVERparams.maxit;  
   x0 = SOLVERparams.x0; 
order = SOLVERparams.order;
    A = SOLVERparams.A;
    b = SOLVERparams.b;
  x_d = SOLUTIONparams.x_d;
% x_ex = SOLUTIONparams.x_ex;
  Nm1 = MESHparams.N-1; 
   M2 = MESHparams.M/2; 

 
%% Construct iteration matrix of the multiplicative Schwarz method

% Notation conflict: In the physical problem we have that each local 
%   subdomain consists of M/2 blocks of size (N-1)x(N-1). To be consistent 
%   with the notatoion of the restriction operators in the manuscript, we 
%   now create local variables where (N-1) is referred to as N and (M/2) 
%   is referred to as (m+1): 
    
   N = Nm1;     n = M2;       m = n-1; 

% We construct the restriction operators with the notation of manuscritpt:  
   I  = eye(N*(m+1)); 
   Z  = zeros(N*(m+1), N*m);

% Restriction operators 
  R1  = [I Z];   
  R2  = [Z I];   

% Approximate inverses
  P1  = R1'*((R1*A*R1')\R1)*A;      
  P2  = R2'*((R2*A*R2')\R2)*A;     

  P1x = R1'*((R1*A*R1')\R1)*b; 
  P2x = R2'*((R2*A*R2')\R2)*b;

% Projections
  Q1 = eye(size(A))-P1; 
  Q2 = eye(size(A))-P2; 

% Construction of iteration matrix with the order of application of
% projections given by the user
if (order==1)  
    T  = Q2*Q1;  
    c = P1x + P2x - P2*P1x;
    
else if (order==2)
    T = Q1*Q2;
    c = P1x + P2x - P1*P2x;
    
    end
    
end 

%% Apply the multiplicative Schwarz method and compute error using

        x_s  = x0;
err_s_inf(1) = norm( x_d - x_s, 'inf');
   err_s_2(1)= norm( x_d - x_s, 2);
   
% err_s_inf(1) = norm( x_ex - x_s, 'inf');
%    err_s_2(1)= norm( x_ex - x_s, 2);

for i = 1:maxit
    res_s(i) = norm(b - A*x_s);
    
         x_s = T*x_s + c;
    
    err_s_inf(i+1) = norm(x_d - x_s, 'inf');
      err_s_2(i+1) = norm(x_d - x_s, 2);
%     err_s_inf(i+1) = norm(x_ex - x_s, 'inf');
%       err_s_2(i+1) = norm(x_ex - x_s, 2);
      
    if err_s_inf(i+1) <= tol, break; end; 
    
end


%%  Store resulting values in SOLUTIONparams data structure

    SOLUTIONparams.x_s   = x_s; 
    SOLUTIONparams.res_s = res_s;
    SOLUTIONparams.err_s_inf = err_s_inf;
    SOLUTIONparams.err_s_2 = err_s_2; 
    SOLUTIONparams.c     = c; 
    SOLUTIONparams.T     = T; 
    SOLUTIONparams.Q1    = Q1; 
    SOLUTIONparams.Q2    = Q2;
    SOLUTIONparams.R1    = R1;
    SOLUTIONparams.R2    = R2;

end