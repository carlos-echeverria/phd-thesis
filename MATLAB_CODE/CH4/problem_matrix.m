function [A] = problem_matrix(Achoice, M, N)

switch Achoice
    
    case 1 %example 4.15a
        
              A = [ 4  -2  -1   1;
                   -2   4   0  -1;
                   -1   0   4  -2;
                    1  -1  -2   4];

    case 2%example 4.15b
        
              A = [ 4.0 -2.0 -0.5  0.5;
                  -2.0  5.0 -1.4 -0.5;
                  -0.5  0.0  4.0 -2.0;
                   0.5 -0.5 -2.0  4.0 ];
               
    case num2cell(3:10)%examples 4.16a,4.16b,4.17,4.18
    
     % Definition of the EQparmas data structure and its subfields
   
    if (Achoice == 3||Achoice == 4||Achoice == 5)         
         EQ.epsi  = 1e-0;
    elseif (Achoice == 6)
         EQ.epsi  = 1e-4;
    else
        EQ.epsi = input('Please specify the value of the perturbation parameter. epsi = ');
    end
   
     EQ.wx  = 0;
     EQ.wy  = 1;
   EQ.beta  = 0;
     EQ.ph1 = 0;
     EQ.ps1 = 0;
     EQ.ph2 = 1;
     EQ.ps2 = 1;

   
%% Definition of the MESHparmas data structure and its subfields

    MESH.N  = N;  %  number of intervals of the grid in x-dir (must be eve1n)
    MESH.M  = M;  %  number of intervals of the grid in y-dir (must be even)
     MESH.scaled = 1;
     MESH.disc   = 1;

     if mod(MESH.N, 2) ~= 0;
         error('The number of intervals in the x-direction is not even.');
     end

      if mod(MESH.M, 2) ~= 0;
         error('The number of intervals in the y-direction is not even.');
     end

%% Definition of the SOLVERparmas data structure and its fields

       K =(MESH.N-1)*(MESH.M-1); % which one?

       SOLVER.tol   = 1e-12;
       SOLVER.maxit = K;
       SOLVER.x0    = zeros(K,1);
       SOLVER.order = 1; % order of application of the projections in the

%% Read off initial prameters from input data structure

        epsi = EQ.epsi;
          wx = EQ.wx;
          wy = EQ.wy;
           N = MESH.N;
           M = MESH.M;


      %% Construct mesh

      % x-direction:
          tau_x = min(0.5,2*epsi*log(N)/wx);    %tau_x = 0.5;
            xi  = zeros(N+1,1);

          for i=0:N/2,
              xi(i+1) = 2*i*(1.0-tau_x)/N;
          end

          for i=N/2+1:N,
              xi(i+1) = 1.0 - 2*(N-i)*tau_x/N;
          end

          Hx = 2*(1-tau_x)/N;
          hx = 2*tau_x/N;

      % y-direction:
          tau_y = min(0.5,2*epsi*log(M)/wy);    %tau_y = 0.5;
            yi  = zeros(M+1,1);

          for i=0:M/2,
              yi(i+1) = 2*i*(1.0-tau_y)/M;
          end

          for i=M/2+1:M,
              yi(i+1) = 1.0 - 2*(M-i)*tau_y/M;
          end

          Hy = 2*(1-tau_y)/M;
          hy = 2*tau_y/M;


      %% Save prameters of constructed mesh to update initial data structure

          MESH.xi = xi;
          MESH.yi = yi;
          MESH.Hx = Hx;
          MESH.hx = hx;
          MESH.Hy = Hy;
          MESH.hy = hy;
          MESH.tau_x = tau_x;
          MESH.tau_y = tau_y;
          
%           fprintf('new problem defined.\n')
%           fprintf('constructing matrix...')

         [SOLVER] = Shishkin2D_Kron_get_A(SOLVER,MESH,EQ);

           A=full(SOLVER.A);
%            figure(77),spy(A)

end

end
