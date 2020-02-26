function [SOLUTIONparams] = block_compute_T_c_inexact(SOLVERparams, MESHparams, SOLUTIONparams)
% BLOCK_COMPUTE_T_C_INEXACT is a function that computes the iteration
% matrix as well as the correction vector (together they build the pre-
% conditioned system) of the mutiplicative Schwarz method with only two
% local subdomain problems.
%
%   function call: 
%
% [SOLUTIONparams] = block_compute_T_c_inexact(MESHparams,SOLVERparams,SOLUTIONparams)
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
% Edited  by C.E. on September 24, 2019.
% Edited  by C.E. on February  25, 2020.


%% Extract needed values from data structures

    A = SOLVERparams.A;
    b = SOLVERparams.b;
   
     order = SOLVERparams.order;
   inexact = SOLVERparams.inexact;
 local_tol = SOLVERparams.local_tol;

  Nm1 = MESHparams.N-1;
  Mm1 = MESHparams.M-1;
   M2 = MESHparams.M/2;


%% Construct iteration matrix of the multiplicative Schwarz method

% Notation conflict: In the physical problem we have that each local
%   subdomain consists of M/2 blocks of size (N-1)x(N-1). To be consistent
%   with the notatoion of the restriction operators in the manuscript, we
%   now create local variables where (N-1) is referred to as N and (M/2)
%   is referred to as (m+1):

   N = Nm1;   M = Mm1;  n = M2;       m = n-1;

   
% We construct the restriction operators with the notation of manuscritpt:
   I  = eye(N*(m+1));
   Z  = zeros(N*(m+1), N*m);

   
% Define restriction operators and local subdomain problems
  R1  = [I Z];
  R2  = [Z I];

   A1 = R1*A*R1';
   A2 = R2*A*R2';
   
   
% Compute exact projections for RHS consistency vector     
  P1  = R1'*(A1\R1)*A;     % approx. inverse 1st domain (exact)
  P2  = R2'*(A2\R2)*A;     % approx. inverse 2nd domain (exact)

  P1b = R1'*(A1\R1)*b;     % projections with P1 (exact)
  P2b = R2'*(A2\R2)*b;     % projections with P2 (exact)
      

% Obtain B and C blocks from transition layer:  
   [ blkA ] = partition_block_matrix(A,N,M);   

       Imp1 = eye(n);
         e1 = Imp1(:,1);
       emp1 = Imp1(:,n);

          B = blkA{n,n+1};
         BB = kron(emp1,B);

          C = blkA{n,n-1};
         CC = kron(e1,C);      

%% Construct iteration matrix T
  
  if (inexact == 0) % exact local solves
      
           
      A1invB = A1\BB;
       Q1_ex = [zeros(N*(m+1),N*(m+1)), -A1invB, zeros(N*(m+1),N*(m-1));
                zeros(N*m,N*(m+1)), eye(N*m)];
            
          Q1 = eye(size(A))-P1;
               
      A2invC = A2\CC;
       Q2_ex = [eye(N*m), zeros(N*m,N*(m+1));
                zeros(N*(m+1),N*(m-1)),-A2invC,zeros(N*(m+1),N*(m+1))];
            
          Q2 = eye(size(A))-P2;
          
        % Construction of iteration matrix with the order of application of
        % projections given by the user
        if (order==1)
%           T = Q2*Q1;
            T = Q2_ex*Q1_ex;
            c = P1b + P2b - P2*P1b;

        else if (order==2)
%           T = Q1*Q2;
            T = Q1_ex*Q2_ex;
            c = P1b + P2b - P1*P2b;

            end

        end

        I_N = eye(size(A));
          P = I_N - T;
        
    SOLUTIONparams.c     = c;
    SOLUTIONparams.T     = T;
    SOLUTIONparams.P     = P;

  else if (inexact==1) % inexact local solves using GMRES
    
      % Use the action of the matrices Q1 and Q2 on a vector:
      
      A1invf = @(x) gmres(A1,BB*x(N*m+N+1:N*m+N+N), [], local_tol, size(A1,1)-1);
%       A1invf = @(x) A1\BB*x(N*m+N+1:N*m+N+N);
      Q1_fun = @(x)[-A1invf(x); x(N*(m+1)+1:N*(2*m+1))];

      A2invf = @(x) gmres(A2,CC*x(N*m-N+1:N*m), [], local_tol, size(A2,1)-1);
%       A2invf = @(x) A2\CC*x(N*m-N+1:N*m);
      Q2_fun = @(x)[x(1:N*m); -A2invf(x)];


        % Construction of iteration matrix:
        if (order==1)
            T_fun = @(x) Q2_fun(Q1_fun(x));  % iteration matrix T = Q2Q1
                c = P1b + P2b - P2*P1b;      % exact

        else if (order==2)
            T_fun = @(x) Q1_fun(Q2_fun(x));  % iteration matrix T = Q1Q2
                c = P1b + P2b - P1*P2b;      % exact
            
            end

        end
        
     P_fun = @(x) x-T_fun(x);  % I-T    
   
     SOLUTIONparams.c     = c;
     SOLUTIONparams.T     = T_fun;
     SOLUTIONparams.P     = P_fun;
  
      end
     
  end


end
