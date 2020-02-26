function [ x_s, res_s, err_s, c, T, P1, P2, R1, R2 ] = Multiplicative_Schwarz( A,b,maxit,tol,order,x0,x,inexact, local_tol)
% Mutiplicative Schwarz method for solving the Shishkin mesh discretized 
% problem.
%
%   function call: 
%
%   [ x_s, res_s, err_s, c, T, P1, P2, R1, R2 ] = Multiplicative_Schwarz(A,b,maxit,tol,order,x0,x)
%
%   input:
%           A: matrix from discretization.
%           b: right hand side of linear system.
%       maxit: maximum number of iterations for the method.
%         tol: convergence tolerance (on error-infinity norm) for the method.
%       order: restriction operator to first domain.
%          x0: initial approximation.
%           x: "exact" solution x = A\b;
%
%   output:
%           x_s: computed approximate solution of linear system. 
%         res_s: vector with the values of the residulas at each iteration.
%         err_s: infinity norms of the errors at each iteration.
%             c: correction vector.
%             T: Multiplicative Scwarz iteration matrix.
%            Pi: Projection  operator to the i-th domain.
%            Ri: Restriction operator to the i-th domain.
%
% Written by Carlos Echeverr?a on January 31, 2016.
% Edited by C.E. on Sept 8, 2016.
% Edited by Joerg Liesen on September 12, 2016.


 N  = length(A)+1;
 n  = N/2;
I_n = eye(n); 
 Z  = zeros(n,n-1);

R1  = [I_n Z];                   % Restriction operator to 1st domain
R2  = [Z I_n];                   % Restriction operator to 2nd domain

if (inexact == 0) 
    
    %% Calculate with exact solves of local problems
    
    P1  = R1'*((R1*A*R1')\R1)*A;     % approx. inverse 1st domain (exact)
    P2  = R2'*((R2*A*R2')\R2)*A;     % approx. inverse 2nd domain (exact)

    P1x = R1'*((R1*A*R1')\R1)*b;     % projections with P1 (exact)
    P2x = R2'*((R2*A*R2')\R2)*b;     % projections with P2 (exact)

else if (inexact==1)
    %% Calculate with inexact solves of local problems with GMRES

    Ah = R1*A*R1';
    AH = R2*A*R2';

    Ih = eye(size(Ah));
    IH = eye(size(AH));

 maxit_g = size(Ah,1);
    x_0 = zeros(size(Ah,1),1);
  Ahinv = zeros(size(Ah));
  AHinv = zeros(size(AH));

    for i = 1:length(Ih)

                       eih = Ih(:,i);
                       eiH = IH(:,i);
    [Ahinv_i,~,~,~,res_Ah] = gmres(Ah,eih,[],local_tol,maxit_g,[],[],x_0);
    [AHinv_i,~,~,~,res_AH] = gmres(AH,eiH,[],local_tol,maxit_g,[],[],x_0);
                Ahinv(:,i) =  Ahinv_i; 
                AHinv(:,i) =  AHinv_i; 
    end

%     norm(Ah*Ahinv-eye(size(Ah)))
%     norm(AH*AHinv-eye(size(AH)))

% s_Ah=size(Ah)
% s_Ahinv=size(Ahinv)
% s_A=size(A)
% s_R1T=size(R1')

    P1  = R1'*(Ahinv)*R1*A;     % approx. inverse 1st domain (inexact)
    P2  = R2'*(AHinv)*R2*A;     % approx. inverse 2nd domain (inexact)

    P1x = R1'*(Ahinv)*R1*b;     % projections with P1 (inexact)
    P2x = R2'*(AHinv)*R1*b;     % projections with P2 (inexact)
    end
end


%% Multiplicative Schwarz method

Q1 = eye(size(A))-P1;    % comlementary projections
Q2 = eye(size(A))-P2;    

if(order==1) 
    T  = Q2*Q1;  % iteration matrix
    c = P1x + P2x - P2*P1x;
else if (order==2)
    T = Q1*Q2;
    c = P1x + P2x - P1*P2x;
    end
end 
 
x_s     = x0;
err_s(1)= norm( x - x_s, 'inf');

for i=1:maxit
    res_s(i) = norm(b - A*x_s);
    x_s = T*x_s + c;
    err_s(i+1) = norm(x - x_s, 'inf');
    if err_s(i+1) <= tol, break; end;       
end
x_s=[0,x_s',0]; 

end

