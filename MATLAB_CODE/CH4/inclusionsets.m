function [ ] = inclusionsets(Achoice)
%INCLUSIONSETS generates the eigenvalue inclusion sets of the given by the
% contour lines FV(z) = ELN(z) = 1, where:
%
%         FV(z) = ||A_{1,2}||*||(A_{1,1}-zI)^{-1}||
%        ELN(z) = ||(A_{1,1}-zI)^{-1}*A_{1,2}||
%
% for the specific cases where A is a 4 x 4 matrix partitioned into blocks
% of size 2x2, given by:
%
%           [ 4  -2 | -1   0]                [   4    -2  | -0.5   0.5]
%       A = [-2   4 |  0  -1]      or    A = [  -2     5  | -1.4  -0.5]
%           [-------|-------]                [------------|-----------]
%           [-1   0 |  4  -2]                [ -0.5    0  |   4    -2 ]
%           [ 0  -1 | -2   4]                [  0.5  -0.5 |   2     4 ]
%
% Written by Carlos Echeverria on April 30, 2018. Last edited on 24.02.20.


%% create a block diagonally dominant matrix A

% Achoice = matA;   % choose matrix to be analyzed

if (Achoice == 1||Achoice == 2)
    M =2;
    N =2;
    I = eye(M);

elseif (Achoice == 3||Achoice == 4||Achoice == 5)
    M=4;
    N=6;
    I = eye(M);
    
elseif (Achoice == 6)
    M = 8;      % choose number of blocks
    N = 16;       % choose block size
else
    M = input('Please specify the number of blocks. M = ');
    N = input('Please specify the size of the blocks. N = ');
    
    M=M+1;
    N=N+1;
end

fprintf('new problem defined.\nconstructing matrix...')

[A] = problem_matrix(Achoice, M, N);

% figure(11), clf,
%      spy(A, 20),hold on, grid on,
%      title('Sparsity pattern of the matrix A ','FontSize',16),
%      set(gca,'FontSize',16),

% InclusionSetMatrix = full(A);
fprintf('done.\npartitioning matrix...')

%% partition matrix 'A' into NxN blocks of size MxM and store them in cells

if (Achoice == 3||Achoice == 4||Achoice == 5||Achoice == 6)
    M =M-1;
    N =N-1;
    I = eye(N);
end


[ blkA ] = PartitionMatrixMN( A, N, M);

if(Achoice == 4)    
    blk12=blkA{1,2};
    blk23=blkA{2,3};
    blk12(1,5)=1e+1;
    blk23(1,5)=1e+1;
    blkA{1,2}=blk12;
    blkA{2,3}=blk23;
    
    
    blk21=blkA{2,1};
    blk32=blkA{3,2};
    blk21(5,1)=1e+1;
    blk32(5,1)=1e+1;
    blkA{2,1}=blk21;
    blkA{3,2}=blk32;
    
elseif (Achoice == 5)
    %makes off-diag blocks random diagonal
    blkA{1,2}=diag(ceil(10*rand(5,1)));
    blkA{2,3}=diag(ceil(10*rand(5,1)));
    
    blkA{2,1}=diag(ceil(10*rand(5,1)));
    blkA{3,2}=diag(ceil(10*rand(5,1)));
    
    %makes off-diag blocks tridiag
%    blkA{1,2}=full(gallery('tridiag',N,-1.6,-16,60));
%    blkA{2,3}=full(gallery('tridiag',N,-1.6,-16,+60));
%    
%    blkA{2,1}=full(gallery('tridiag',N,-0.1,-20,+10));
%    blkA{3,2}=full(gallery('tridiag',N,-0.1,-20,+10));

end
    
    A = cell2mat(blkA);
    
     figure(11), clf,
      spy(A),
      title('Spsrsity pattern of the matrix A ','FontSize',16),
      set(gca,'FontSize',16),



fprintf('done.\n')
if (Achoice==1||Achoice==2||Achoice==3||Achoice==4||Achoice==5)
    fprintf('The matrix to be handled is:\n')

         A = full(A)
    [~, D] = eig(A); % compute eigenvalues of chosen matrix
     d = diag(D);
%    e = zeros(size(d));
      evals = sort(d)
else
             A = full(A);
        [~, D] = eig(A); % compute eigenvalues of chosen matrix
             d = diag(D);
         evals = sort(d)
%    e = zeros(size(d));
end

     
fprintf('evaluating submultiplicativity condition of block rows...')

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

fprintf('done.\n')

if (ELN>=FV)
    where=find(ELN>=FV);
    warning('Matrix A does not fulfill strict submultiplicativity condition on block-row %d.\n', where);
    warning('Inclusion sets from ELN might not be contained in the sets FV.')
end

%% define BDiDo functions FV and ELN

fprintf('recursively defining complex Block Diag-Dom. functions of FV and ELN...')

%fELN{k} will be the k-th Gershgorin set.
for k=1:M
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

fprintf('done.\ncreating computational mesh...')

%% mesh the complex plane and evaluate BDiDo functions

if (Achoice==1||Achoice==2)
    
    xmin = -2.0;
    xmax = 10;
    ymin = -2.2;
    ymax =  2.2;
    
    Nptsx = 50; % NUmber of points between xmin and xmax
    Nptsy = 50; % NUmber of points between ymin and ymax
    
elseif (Achoice==3||Achoice==4)
    
     xmin = -20;
     xmax = 220;
     ymin = -50;
     ymax = 50;
    
    Nptsx = 60; 
    Nptsy = 60; 

elseif (Achoice==5)
    
     xmin = 0;
     xmax = 200;
     ymin = -25;
     ymax =  25;
    
    Nptsx = 50; 
    Nptsy = 50; 
    
elseif (Achoice==6)
    
     xmin = -5e3;
     xmax =  6e4;
     ymin = -3e4;
     ymax =  3e4;
    
    Nptsx = 20; 
    Nptsy = 20;
    
else
     xmin = min(real(d))-2.5*abs(max(real(d)))/2;
     xmax = max(real(d))+0.9*abs(max(real(d)))/2;
     ymin = min(imag(d))-1.1*abs(max(real(d)))/2;
     ymax = max(imag(d))+1.1*abs(max(real(d)))/2;
     
     Nptsx = 20; 
     Nptsy = 20; 
end

    
    x = linspace(xmin, xmax, Nptsx);
    y = linspace(ymin, ymax, Nptsy);
[X,Y] = meshgrid(x,y);

fprintf('done.\nevaluating functions on computational mesh of size %dx%d...', Nptsx, Nptsy)

%% evaluate function on discrete complex plane
for k=1:M
    for i=1:size(x,2)
        for j=1:size(y,2)
            EELN{k}(j,i)= fELN{k}(x(i),y(j));

            FFV{k}(j,i)= fFV{k}(x(i),y(j));
        end
    end
end
fprintf('done.\nplotting results...');

%% plot complex functions and eigenvalue inclusion regions

if (Achoice==1||Achoice==2)
    plotSets_ex1( X, Y, FFV, EELN, d);
    fprintf('done.\nend.\n');
else
    plotSets( X, Y, FFV, EELN, d);
    fprintf('done.\nend.\n');
          matrix_size = size(A) 
      number_of_evals =size(d,1)
end


end
