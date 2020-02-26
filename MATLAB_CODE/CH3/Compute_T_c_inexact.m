function [ c, T, B, P1_ex, P2_ex, R1, R2 ] = Compute_T_c_inexact( A,b,order, inexact, local_tol)
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
% Written by Carlos Echeverria on January 31, 2016.
% Edited  by C.E. on Sept 8, 2016.
% Edited  by Joerg Liesen on September 12, 2016.
% Edited  by Daniel B. Szyld from the code: Multiplicative_Schwarz.m on December 30, 2016.
% Edited  by C.E. on October 16, 2019
% Edited  by C.E. on February 20, 2020.


    N  = length(A)+1;
    n  = N/2;
   I_n = eye(n); 
    Z  = zeros(n,n-1);

    R1 = [I_n Z];           % Restriction operator to 1st domain
    R2 = [Z I_n];           % Restriction operator to 2nd domain


    A1 = R1*A*R1';
    A2 = R2*A*R2';
    
%% Calculate local problems with exact solves

P1_ex  = R1'*((R1*A*R1')\R1)*A;     % approx. inverse 1st domain (exact)
P2_ex  = R2'*((R2*A*R2')\R2)*A;     % approx. inverse 2nd domain (exact)

P1b_ex = R1'*((R1*A*R1')\R1)*b;     % projections with P1 (exact)
P2b_ex = R2'*((R2*A*R2')\R2)*b;     % projections with P2 (exact)

%% Multiplicative Schwarz method

Q1_ex1 = eye(size(A))-P1_ex;    % comlementary projections
Q2_ex1 = eye(size(A))-P2_ex;    

I1 = eye(size(A1));

b_entry = A(n,n-1);
en = I1(:,n);
A1inv_ex = b_entry*A1\en;

Q1_ex = [ zeros(n,n), -A1inv_ex, zeros(n,n-2);
           zeros(n-1,n), eye(n-1)];


c_entry = A(n,n+1);
    e1 = I1(:,1);
A2inv_ex = c_entry*A2\e1;
Q2_ex = [eye(n-1), zeros(n-1,n);
          zeros(n,n-2),-A2inv_ex, zeros(n,n)];

norm(Q1_ex1-Q1_ex)
norm(Q2_ex1-Q2_ex)


if(order==1) 
    T_ex  = Q2_ex*Q1_ex;  % iteration matrix
    c_ex = P1b_ex + P2b_ex - P2_ex*P1b_ex;
else if (order==2)
    T_ex = Q1_ex*Q2_ex;
    c_ex = P1b_ex + P2b_ex - P1_ex*P2b_ex;
    end
end 


B_ex = eye(size(A)) - T_ex;

T = T_ex;
B = B_ex;
c = c_ex;

if (inexact == 1) 
    %% Calculate with inexact solves of local problems with GMRES
    

    I1 = eye(size(A1));
    
       b_entry = A(n,n-1);
            en = I1(:,n);
    A1inv_inex = gmres(A1,b_entry*en, [], local_tol, size(A1,1)-1);
       Q1_inex = [ zeros(n,n), -A1inv_inex, zeros(n,n-2);
                   zeros(n-1,n), eye(n-1)];
 
    
       c_entry = A(n,n+1);
            e1 = I1(:,1);
    A2inv_inex = gmres(A2,c_entry*e1, [], local_tol, size(A2,1)-1);
       Q2_inex = [eye(n-1), zeros(n-1,n);
                  zeros(n,n-2),-A2inv_inex, zeros(n,n)];
  
              
    if(order==1) 
        T_inex = Q2_inex*Q1_inex;  % iteration matrix
    else if (order==2)
        T_inex = Q1_inex*Q2_inex;
        end
    end 
    
    norm(T_ex-T_inex)
%       b = A(n+1,n);
%      e1 = I1(:,1);
%     A1invf = gmres(A1,b*e1, [], local_tol, size(A1,1)-1);
%     Q1_fun = @(x) [-A1invf; x(n+1:N-1)];
% 
%     
%      c = A(n+1,n+2);
%     en = I1(:,end);
%     A2invf = gmres(A2,c*en, [], local_tol, size(A2,1)-1);
%     Q2_fun = @(x) [x(1:n); -A2invf];
% 
% 
%     if(order==1) 
%         T_fun  = @(x) Q2_fun(Q1_fun(x));  % iteration matrix
%     else if (order==2)
%         T_fun = @(x) Q1_fun(Q2_fun(x));
%         end
%     end 


    B_inex = eye(size(A,1)) - T_inex;

    T = T_inex;
    B = B_inex;
    
end



%% Calculate error
% 
% x0 = zeros(N-1,1);
% x_s     = x0;
% err_s(1)= norm( x - x_s, 'inf');
% 
% for i=1:maxit
%     res_s(i) = norm(b - A*x_s);
%     x_s = T*x_s + c;
%     err_s(i+1) = norm(x - x_s, 'inf');
%     if err_s(i+1) <= tol, break; end;       
% end
% x_s=[0,x_s',0]; 



end


















% function [ c, T, B, P1, P2, R1, R2 ] = Compute_T_c_inexact( A,b,order, inexact, local_tol)
% % Computation of preconditioned matrix with 
% % Mutiplicative Schwarz method for solving the Shishkin mesh discretized 
% % problem.
% %
% %   function call: 
% %
% %   [ c, T, P1, P2, R1, R2 ] = Compute_T_c(A,b,order)
% %
% %   input:
% %           A: matrix from discretization.
% %           b: right hand side of linear system.
% %       order: restriction operator to first domain.
% %
% %   output:
% %             c: correction vector.
% %             T: Multiplicative Scwarz iteration matrix.
% %             B: I-T
% %            Pi: Projection  operator to the i-th domain.
% %            Ri: Restriction operator to the i-th domain.
% %
% % modified by DBS From the code: Multiplicative_Schwarz.m on 30 Dec 2016
% %
% % Written by Carlos Echeverr?a on January 31, 2016.
% % Edited by C.E. on Sept 8, 2016.
% % Edited by Joerg Liesen on September 12, 2016.
% 
% 
%  N  = length(A)+1;
%  n  = N/2;
% I_n = eye(n); 
%  Z  = zeros(n,n-1);
% 
% R1  = [I_n Z];                   % Restriction operator to 1st domain
% R2  = [Z I_n];                   % Restriction operator to 2nd domain
% 
% 
% if (inexact == 0) 
%     
%     %% Calculate with exact solves of local problems
%     
%     P1  = R1'*((R1*A*R1')\R1)*A;     % approx. inverse 1st domain (exact)
%     P2  = R2'*((R2*A*R2')\R2)*A;     % approx. inverse 2nd domain (exact)
% 
%     P1x = R1'*((R1*A*R1')\R1)*b;     % projections with P1 (exact)
%     P2x = R2'*((R2*A*R2')\R2)*b;     % projections with P2 (exact)
% 
% else if (inexact==1)
%     %% Calculate with inexact solves of local problems with GMRES
% 
%     Ah = R1*A*R1';
%     AH = R2*A*R2';
% 
%     Ih = eye(size(Ah));
%     IH = eye(size(AH));
% 
%  maxit_g = size(Ah,1);
%     x_0 = zeros(size(Ah,1),1);
%   Ahinv = zeros(size(Ah));
%   AHinv = zeros(size(AH));
% 
%     for i = 1:length(Ih)
% 
%                        eih = Ih(:,i);
%                        eiH = IH(:,i);
%     [Ahinv_i,~,~,~,res_Ah] = gmres(Ah,eih,[],local_tol,maxit_g,[],[],x_0);
%     [AHinv_i,~,~,~,res_AH] = gmres(AH,eiH,[],local_tol,maxit_g,[],[],x_0);
%                 Ahinv(:,i) =  Ahinv_i; 
%                 AHinv(:,i) =  AHinv_i; 
%     end
% 
% 
%     P1  = R1'*(Ahinv)*R1*A;     % approx. inverse 1st domain (inexact)
%     P2  = R2'*(AHinv)*R2*A;     % approx. inverse 2nd domain (inexact)
% 
%     P1x = R1'*(Ahinv)*R1*b;     % projections with P1 (inexact)
%     P2x = R2'*(AHinv)*R1*b;     % projections with P2 (inexact)
%     end
% end
% 
% 
% %% Multiplicative Schwarz method
% 
% Q1 = eye(size(A))-P1;    % comlementary projections
% Q2 = eye(size(A))-P2;    
% 
% if(order==1) 
%     T  = Q2*Q1;  % iteration matrix
%     c = P1x + P2x - P2*P1x;
% else if (order==2)
%     T = Q1*Q2;
%     c = P1x + P2x - P1*P2x;
%     end
% end 
% 
% 
% B = eye(size(A)) - T;
%  
% 
% %% Calculate error
% % 
% % x0 = zeros(N-1,1);
% % x_s     = x0;
% % err_s(1)= norm( x - x_s, 'inf');
% % 
% % for i=1:maxit
% %     res_s(i) = norm(b - A*x_s);
% %     x_s = T*x_s + c;
% %     err_s(i+1) = norm(x - x_s, 'inf');
% %     if err_s(i+1) <= tol, break; end;       
% % end
% % x_s=[0,x_s',0]; 
% 
% 
% 
% end