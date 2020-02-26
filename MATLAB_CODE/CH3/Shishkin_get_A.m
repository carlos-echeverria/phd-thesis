function [ A,D ] = Shishkin_get_A( eps,alpha,beta,N,disc,scaled )
%Shishkin_get_A Constructs a matrix corresponding to the discretization of
% the 1D constant-coefficient convection-diffusion equation on a Shishkin 
% mesh.
%
%   This function generates a finite difference discretization of the 
%   second order elliptic partial differential equation:
%
%             -eps*u'' + alpha*u' + beta*u = f, u(0)=0, u(1)=0,
%
%   on the interval [0,1] discretized by a Shishkin mesh with N intervals.
%   The discretization can be chosen between an upwind or a central 
%   approximation of the first derivatives. A scaling factor can also be
%   included that correctly scales the matrix of eigenvectors.
%
%   function call:
%
%       [ A,D ] = Shishkin_get_A( eps,alpha,beta,N,disc,scaled )
%
%   input: 
%             eps: scalar diffusion  coefficient.  
%           alpha: scalar convection coefficient.
%            beta: scalar reaction   coefficient.
%               N: the number of intervals conformng the grid.
%            disc: specifies the type of discretization.
%                  can be 1 for upwind or 2 for central differencing.
%          scaled: can be 1 for a scaled matrix or 0 for unscaled.
%
%   output:
%               A: System matrix.
%               D: Scaling matrix.
%
% Written by Carlos Echeverria on January 31, 2016.
% Edited by C.E. on September 8, 2016.
% Edited by Joerg Liesen on September 12, 2016.
% Edited by C.E. on February 20, 2020.

format long
n = N/2;
tau = min(0.5, 2*eps*log(N)/alpha);
H = 2*(1-tau)/N; h = 2*tau/N;

%% Choose discretization type

switch disc
    case 1 % simple upwind discretization
        v_W = -(eps/(H^2)) - (alpha/H); 
        v_C = ((2*eps)/(H^2)) + (alpha/H) + beta;
        v_E =   -eps/(H^2);

        z_W = ((-2*eps)/((H^2)+(h*H)))-(alpha/H); 
        z_C = ((2*eps)/(h*H))+(alpha/H)+beta;
        z_E = -2*eps /(h*H+h^2);

        w_W = -(eps/(h^2))    - (alpha/h); 
        w_C = ((2*eps)/(h^2)) + (alpha/h)+beta;
        w_E = -eps/(h^2);
        
        if(scaled==1) 
            D = diag([(H/alpha)*ones(n-1,1);...
                      (h*H)/(2*eps);((h*h)/eps)*ones(n-1,1)]);
        else
            D = eye(N-1);
        end
    
    case 2 % central diference discretzation
        v_W = -(eps/(H^2))-(alpha/(2*H)); 
        v_C = ((2*eps)/H^2)+beta;
        v_E = -(eps/(H^2))+(alpha/(2*H));

        z_W = ((-2*eps)/(H^2+h*H))-(alpha/(H+h)); 
        z_C = ((2*eps)/(h*H))+beta;
        z_E = -2*eps/(h^2+h*H)+(alpha/(h+H));

        w_W = -(eps/(h^2))-(alpha/(2*h)); 
        w_C = ((2*eps)/(h^2))+beta;
        w_E = -eps/(h^2)+alpha/(2*h);
        
        if(scaled==1) 
            D = diag([((2*H)/alpha)*ones(n-1,1);((h*H)+(h*h))/(2*eps);...
                       ((h*h)/eps)*ones(n-1,1)]);
        else
            D = eye(N-1);
        end
end 
 

%% Construct system matrix

A  = zeros(N-1,N-1);
e  = ones(n,1);
D1 = spdiags([v_W*e,v_C*e,v_E*e],-1:1,n-1,n-1); % outside layer submatrix
D2 = spdiags([w_W*e,w_C*e,w_E*e],-1:1,n-1,n-1); % inside layer submatrix

A (1:n-1  , 1:n-1  ) = D1;
A (n+1:N-1, n+1:N-1) = D2;
A (n      , n      ) = z_C; 
A (n      , n-1    ) = z_W;
A (n      , n+1    ) = z_E; 
A (n-1    , n      ) = v_E;
A (n+1    , n      ) = w_W;

A = D*A;

end