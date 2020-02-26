function [ A, Aparams ] = getA( Achoice )
%GETA builds an M^2 x M^2 block tridiagonal Toeplitz matrix with
%  tridiagonal blocks.
%
%   function call:
%
%          [ A ] = getA( Aparams, M )
%
%   input:
%
%      Aparams.gamma: subdiagonal of diagonal block
%       Aparams.beta: superdiagonal of diagonal block
%      Aparams.alpha: diagonal of diagonal block
%         Aparams.s1: superdiagonal of superdiagonal block
%         Aparams.s2: subdiagonal of superdiagonal block
%         Aparams.s3: superdiagonal of subdiagonal block
%         Aparams.s4: subdiagonal of superdiagonal block
%                  M: size of the blocks.
%
%   output:
%
%                  A: block tridiagonal matrix.
%
% Written by Carlos Echeverria. Last edited on Nov. 25, 2017.
M=9;
switch Achoice

    case 1, % 2D Laplacian on regular mesh (symmetric).
             Aparams.M = M;
         Aparams.gamma = -1;
         Aparams.beta  = -1;
         Aparams.alpha = abs(Aparams.gamma)+abs(Aparams.beta);
            Aparams.s1 = 0;
            Aparams.s2 = 0;
            Aparams.s3 = 0;
            Aparams.s4 = 0;
       Aparams.Achoice = 1;

    case 2, % 2D Convection Diffusion on regular mesh (unsymmetric).

             Aparams.M = M;
         Aparams.gamma = -1.100000000000000e+02;
         Aparams.beta  = -9.999999999999999e+01;
         Aparams.alpha =  2.099999999999999e+02;
            Aparams.s1 = 0;
            Aparams.s2 = 0;
            Aparams.s3 = 0;
            Aparams.s4 = 0;
       Aparams.Achoice = 2;

    case 3, % Random block tridiagonal matrix which is BDDom. (diagonal off-diag. blocks)

             Aparams.M = M;
         Aparams.gamma = -1;
         Aparams.beta  = -1;
         Aparams.alpha = abs(Aparams.gamma)+abs(Aparams.beta);
            Aparams.s1 = 0;
            Aparams.s2 = 0;
            Aparams.s3 = 0;
            Aparams.s4 = 0;
       Aparams.Achoice = 3;

    case 4, % Unsymmetric matrix which fulfills submultiplicativity cond. striclty.

             Aparams.M = M;
         Aparams.gamma = -2;
         Aparams.beta  = -2;
         Aparams.alpha = 5;
            Aparams.s1 = .01;
            Aparams.s2 = -1;
            Aparams.s3 = .01;
            Aparams.s4 = -1;
       Aparams.Achoice = 4;

end

        gamma = Aparams.gamma;
         beta = Aparams.beta;
        alpha = Aparams.alpha;
           s1 = Aparams.s1;
           s2 = Aparams.s2;
           s3 = Aparams.s3;
           s4 = Aparams.s4;

%if alpha ~= abs(beta)+abs(gamma), alpha = 0.5*alpha; end

            I = eye(M);
            T = gallery('tridiag',M,gamma,alpha,beta);
            A = kron(I,T) + kron(T,I);
           D1 = diag(s1*ones(1,M^2-(M+1)),M+1);
           D2 = diag(s2*ones(1,M^2-(M-1)),M-1);
           D3 = diag(s3*ones(1,M^2-(M-1)),-(M-1));
           D4 = diag(s4*ones(1,M^2-(M+1)),-(M+1));
            A = A + D1 + D2 + D3 + D4;

% remove elements to make matrix block tridiagonal:

        for i = 1:M-2
            A(i*M,(i+1)*M+1)=0; %remove unwanted elements from D1
            A((i+1)*M+1,i*M)=0; %remove unwanted elements from D4
        end

% without removing these elements the matrix is still block tridigonal:

        for i = 1:M
            A((i-1)*M+1,i*M)=0; %remove unwanted elements from D2
            A(i*M,(i-1)*M+1)=0; %remove unwanted elements from D3
        end

% randomize entries for examples 2.11 and 2.12
        if Aparams.Achoice==3 || Aparams.Achoice==4,

                R1 = diag(ceil(diag(10*rand(M))));
                R2 = diag(ceil(diag(10*rand(M))));
                 A = kron(R1,eye(M))*A;
        end

figure(1), clf,
     spy(A), set(gca,'FontSize',16),
     title('Spsrsity pattern of the matrix A ','FontSize',16), pause(2);

end
