function [ ] = inclusionsets2D(N,M)
%INCLUSIONSETS generates the eigenvalue inclusion sets of the given by the
% contour lines FV_i(z) = ELN_i(z) = 1, where:
%
%     FV_i(z) = SUM_{j=1,j!=i}^{N}(||A_{ij}||*||(A_{ii}-zI)^{-1}||)
%    ELN_i(z) = SUM_{j=1,j!=i}^{N}(||(A_{ii}-zI)^{-1}*A_{ij}||)
%
% for the specific cases where A comes from the 1D discretization of the Shishkin problem
%
% Written by Carlos Echeverria on October 28, 2019.


%% create a block diagonally dominant matrix A


%% Definition of the EQparmas data structure and its subfields
   EQ.epsi  = 1e-0;
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
fprintf('new problem defined.\n')
fprintf('constructing matrix...')

 [SOLVER] = Shishkin2D_Kron_get_A(SOLVER,MESH,EQ);

 A=SOLVER.A;
 figure(77),spy(A)
fprintf('done\npartitioning matrix...')
 figure(11), clf,
      spy(A),
      title('Spsrsity pattern of the matrix A ','FontSize',16),
      set(gca,'FontSize',16),
%       pause(5);

% InclusionSetMatrix = full(A)
% save('matAA.mat', 'InclusionSetMatrix')

%% partition matrix 'A' into NxN blocks ofsize MxM and store them in cells

M = M-1;
N = N-1;
 I = eye(N);
[ blkA ] = PartitionMatrixMN( A, N, M);
fprintf('done\nevaluating matrix submultiplicativity conditions...')

%% modify off-diagonal block entries (if needed set flag=true)

  flag = false(1);
%   flag = true(1);
if(flag==1)
    
    %changes upper right, lower left corners
%     blk12=blkA{1,2};
%     blk23=blkA{2,3};
%     blk12(1,5)=1e+1;
%     blk23(1,5)=1e+1;
%     blkA{1,2}=blk12;
%     blkA{2,3}=blk23;
%     
%     
%     blk21=blkA{2,1};
%     blk32=blkA{3,2};
%     blk21(5,1)=1e+1;
%     blk32(5,1)=1e+1;
%     blkA{2,1}=blk21;
%     blkA{3,2}=blk32;
    
        %makes off-diag blocks tridiag
%     blkA{1,2}=full(gallery('tridiag',N,-1.6,-16,60));
%     blkA{2,3}=full(gallery('tridiag',N,-1.6,-16,+60));
%     
%     blkA{2,1}=full(gallery('tridiag',N,-0.1,-20,+10));
%     blkA{3,2}=full(gallery('tridiag',N,-0.1,-20,+10));

    %makes off-diag blocks random diagonal
    blkA{1,2}=diag(floor(10*rand(5,1)));
    blkA{2,3}=diag(floor(10*rand(5,1)));
    
    blkA{2,1}=diag(floor(10*rand(5,1)));
    blkA{3,2}=diag(floor(10*rand(5,1)));

    
    A = cell2mat(blkA);
     figure(11), clf,
      spy(A),
      title('Spsrsity pattern of the matrix A ','FontSize',16),
      set(gca,'FontSize',16),
end

%% test for strictness of matrix norm submultiplicativity condition in each block row

ELN = zeros(1,M);
 FV = zeros(1,M);
for i = 1:M-1
    for j = i+1:M
         ELN(i) = ELN(i) + norm((blkA{i,i}\I)*blkA{i,j});
          FV(i) =  FV(i) + norm(blkA{i,i}\I)*norm(blkA{i,j});
    end
end

for i = 2:M
    for j = 1:i-1
       ELN(i) = ELN(i) + norm((blkA{i,i}\I)*blkA{i,j});
        FV(i) = FV(i) + norm(blkA{i,i}\I)*norm(blkA{i,j});
    end
end

% ELN = ELN;
% FV = FV;
if (ELN>FV)
    where=find(ELN>FV);
    warning('Matrix A does not fulfill strict submultiplicativity condition on row %d.\n', where);
    warning('Inclusion sets from ELN might not be contained in the sets FV.\n', where);
end
fprintf('done\nrecursively defining BDiDo functions...')

%% compute eigenvalues of chosen matrix


%fELN{k} will be the k-th Gershgorin set.
for k = 1:M
    fELN{k} = @(x,y) 0; % initialize the functions as zero
     fFV{k} = @(x,y) 0; % initialize the functions as zero
end

for ii = 1:M-1
    for jj = ii+1:M
          fELN_j{jj} = @(x,y) norm((blkA{ii,ii}-(x+1i*y).*I)\blkA{ii,jj});
          fELN{ii} = @(x,y) fELN{ii}(x,y) + fELN_j{jj}(x,y);

          fFV_j{jj} = @(x,y) norm((blkA{ii,ii}-(x+1i*y).*I)\I)*norm(blkA{ii,jj});
          fFV{ii} = @(x,y) fFV{ii}(x,y) + fFV_j{jj}(x,y);
    end
end

for ii = 2:M
    for jj = 1:ii-1
    fELN_j{jj} = @(x,y) norm((blkA{ii,ii}-(x+1i*y).*I)\blkA{ii,jj});
       fELN{ii} = @(x,y) fELN{ii}(x,y) + fELN_j{jj}(x,y);

    fFV_j{jj} = @(x,y) norm((blkA{ii,ii}-(x+1i*y).*I)\I)*norm(blkA{ii,jj});
       fFV{ii} = @(x,y) fFV{ii}(x,y) + fFV_j{jj}(x,y);
    end
end

fprintf('done\ncreating computational mesh...')
%% compute eigenvalues of chosen matrix
    
           A = full(A)
      size(A)

     
   [EVEC, D] = eig(A);
           d = diag(D);

      sort(d)
      size(d)
        
        e = EVEC;
        size(e);
%% mesh the complex plane and evaluate BDiDo functions
%      xmin = min(real(d))-2.5*abs(max(real(d)))/2;
%      xmax = max(real(d))+0.9*abs(max(real(d)))/2;
%      ymin = min(imag(d))-1.1*abs(max(real(d)))/2;
%      ymax = max(imag(d))+1.1*abs(max(real(d)))/2;
     
     xmin = -20;
     xmax = 220;
     ymin = -50;
     ymax = 50;

%      xmin = -0;
%      xmax = 200;
%      ymin = -25;
%      ymax = 25;
     
%eps=1
%      xmin = -200;
%      xmax = 1500;
%      ymin = -200;
%      ymax =  200;


%eps=1e-1
%      xmin = -20;
%      xmax = 160;
%      ymin = -50;
%      ymax =  50;

% %eps=1e-2
%      xmin = -50;
%      xmax = 600;
%      ymin = -300;
%      ymax =  300;


% %eps=1e-3
%      xmin = -100;
%      xmax = 6000;
%      ymin = -3000;
%      ymax =  3000;


%eps=1e-4
%      xmin = -5000;
%      xmax = 60000;
%      ymin = -30000;
%      ymax =  30000;

    Nptsx = 20; % Number of points between xmin and xmax
    Nptsy = 20; % Number of points between ymin and ymax

        x = linspace(xmin, xmax, Nptsx);
        y = linspace(ymin, ymax, Nptsy);
    [X,Y] = meshgrid(x,y);

fprintf('done\nevaluating BdiDo functions on mesh of size %dx%d...', Nptsx, Nptsy)
%% evaluate function on discrete complex plane
for k=1:M
    for i=1:size(x,2)
        for j=1:size(y,2)
            EELN{k}(j,i)= fELN{k}(x(i),y(j));

            FFV{k}(j,i)= fFV{k}(x(i),y(j));
        end
    end
end
fprintf('done\nplotting results...')
%% plot complex functions and eigenvalue inclusion regions
figure(8),clf,
purple = [.49 .18 .55];
blue = [.301 .745 .933];

plot(real(d), imag(d), 'k.', 'MarkerSize', 28), hold on,
linS = {'--','-.',':'};
% linS = fliplr(linS);
for ii = 1:length(FFV)
    txt = sprintf('G_{%d}^{FV}', ii);
    [C,h] = contour(X, Y,  FFV{ii},  [1 1], 'LineWidth', 3, 'Color', purple);

    %clabel(C,h,'FontSize',15);
end
linS = {'-',':','-.','--'};
% linS = fliplr(linS);
for ii = 1:length(EELN)
    txt = sprintf('G_{%d}^{ELN}', ii);
    [C,h] = contour(X, Y,  EELN{ii},  [1 1], 'LineWidth', 3,'Color', blue);
    %clabel(C,h,'FontSize',15);
end

axis equal, grid on,
% leg = legend('eigenvalues', 'G_{1}^{FV}', 'G_{2}^{FV}', 'G_{3}^{FV}', 'G_{1}^{new}', 'G_{2}^{new}', 'G_{3}^{new}');
% set(leg,'FontSize',17,'Location','BestOutside','Orientation','horizontal')
xlabel('Re(z)','FontSize',17),
ylabel('Im(z)','FontSize',17),
% title('Eigenvalue Inclusion Regions','FontSize',17)
set(gca,'FontSize',18),
fprintf('done\nend.\n\n')
end
