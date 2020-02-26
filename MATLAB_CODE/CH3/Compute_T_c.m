function [ c, T, B, P1, P2, R1, R2 ] = Compute_T_c(A,b,order)
% Computation of preconditioned matrix with 
% Mutiplicative Schwarz method for solving the Shishkin mesh discretized 
% problem.
%
%   function call: 
%
%   [ c, T, P1, P2, R1, R2 ] = Compute_T_c(A,b,order)
%
%   input:
%           A: matrix from discretization.
%           b: right hand side of linear system.
%       order: restriction operator to first domain.
%
%   output:
%             c: correction vector.
%             T: Multiplicative Scwarz iteration matrix.
%             B: I-T
%            Pi: Projection  operator to the i-th domain.
%            Ri: Restriction operator to the i-th domain.
%
% modified by DBS From the code: Multiplicative_Schwarz.m on 30 Dec 2016
%
% Written by Carlos Echeverria on January 31, 2016.
% Edited by C.E. on Sept 8, 2016.
% Edited by Joerg Liesen on September 12, 2016.



 N  = length(A)+1;
 n  = N/2;
I_n = eye(n); 
 Z  = zeros(n,n-1);

R1  = [I_n Z];                   % Restriction operator to 1st domain
R2  = [Z I_n];                   % Restriction operator to 2nd domain


P1  = R1'*((R1*A*R1')\R1)*A;     % approx. inverse 1st domain
P2  = R2'*((R2*A*R2')\R2)*A;     % approx. inverse 2nd domain 

P1x = R1'*((R1*A*R1')\R1)*b; 
P2x = R2'*((R2*A*R2')\R2)*b;

M1 = eye(size(A))-P1; M2 = eye(size(A))-P2; % projections

if(order==1) 
    T  = M2*M1;  % iteration matrix
    c = P1x + P2x - P2*P1x;
else if (order==2)
    T = M1*M2;
    c = P1x + P2x - P1*P2x;
    end
end 

I_N = eye(size(A));

B = I_N - T;
 


end

