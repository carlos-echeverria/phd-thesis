% Script for performing the numerical experiments concerning the
% theoretical bounds presented in Chapter 4 of the PH.D. thesis:
%
%  [ Echeverria - Iterative solution of discretized convection
%    diffusion problems -Technische Universitaet Berlin - 2020 ]
%
%  In order to perform the experiments, it is needed to hardcode the
%  input into this file. It consist of three parameters which can be edited
%  in the section "chose parameters of the problem" below.
%
% input:
%
%     Achoice: matrix to be analyzed. It accepts values 1 to 4, where:
%                    1 is the matrix of Example 4.11,
%                    2 is the matrix of Example 4.12,
%                    3 is the matrix of Example 4.13,
%                    4 is the matrix of Example 4.14.
%
%   directinv: flag to decide how to calculate the inverse of A.
%              1 for backslash or 0 for sequences UVXY.
%
%   refinement_step: flag to indicate the desired level of refinement of
%                    the bounds.
%
% subordinate functions:
%
%            getA.m
%            PartitionMatrix.m
%            BlockDiagDomCheck.m
%            CommutativityCheck.m
%            inverseUVXY.m
%            PartitionMatrix.m
%            TausAndOmegas.m  % revised version of omegas
%            TausAndOmegas2.m % unrevised version of omegas
%            generateBounds.m
%            plotBounds.m
%
% output:
%
%        * all the figures for the chosen example.
%        * values of the maximum error for the upper and lower bounds
%        * values of the condition number and deviation from identity
%
% Written by Carlos Echeverria on May 4, 2018.
% Editted by Carlos Echeverria on October 31, 2019.
% Editted by Carlos Echeverria on February 24, 2020.


clearvars, close all

%% choose parameters of the problem

          Achoice = 4; % choose example to be analyzed
        directinv = 0; % choose matrix inversion method
  refinement_step = 8; % choose level of bound refinement (at most M-1)

%% construct matrix 'A' from given parameters

[A, Aparams] = getA(Achoice);
        cond = condest(A);
         % A = importdata('A4.mat') % uncoment for importing matrices from file

%% partition matrix 'A' into blocks and store them in cells

[blkA] = PartitionMatrix(A);

%% check for row block diagonal dominance of chosen 'A'

[BDD, BDDFV] = BlockDiagDomCheck( blkA );

%% check for commutativity between off diagonal and diagonal blocks of 'A':

[commuteAB, commuteAC] = CommutativityCheck( blkA );

 %% compute 'Z', the inverse of 'A', via sequences or direct inversion
 M = Aparams.M;
switch directinv

    case 1
            Z = A\eye(M*M);
    case 0
       [ ZZ ] = inverseUVXY( blkA );
            Z = cell2mat(ZZ);
end

 % check to what extent the computed Z fulfills Z*A = I
 dev  = norm(Z*A-eye(M*M));

%% partition matrix 'Z' into blocks and store them in cells

[blkZ] = PartitionMatrix(Z);

%% build taus und omegas

  [ tauu, omegaa ] = TausAndOmegas( blkA );  % revised version of omegas
[ tauu2, omegaa2 ] = TausAndOmegas2( blkA ); % unrevised version of omegas

%% Generate bounds for the norms of the blocks of 'Z' for all t and calculate errors

 MaxErrU = zeros(1,M-1);
 MaxErrL = zeros(1,M-1);

 for t = 1:refinement_step;

 % generate bounds with revised omegas:
[UBound,LBound,upbound,lowbound,entries,diagentries] = ...
                              generateBounds(blkZ, blkA, tauu, omegaa, t );

            UError = (UBound-entries)./UBound;
            LError = (diagentries-lowbound)./diagentries;

        MaxErrU(t) = max(max(UError));
        MaxErrL(t) = max(LError);

% generate bounds with unrevised omegas:
[UBound2,LBound2,upbound2,lowbound2,entries2,diagentries2] = ...
                            generateBounds(blkZ, blkA, tauu2, omegaa2, t );

            UError2 = (UBound2-entries2)./UBound2;
            LError2 = (diagentries2-lowbound2)./diagentries2;

        MaxErrU2(t) = max(max(UError2));
        MaxErrL2(t) = max(LError2);

        [ ~ ] = plotBounds( UBound, LBound, entries, diagentries, upbound, lowbound, Aparams, MaxErrU, t );


 end

%% output results to the console
format shorte

        condition_number = cond
 deviation_from_identity = dev
           max_err_upper = MaxErrU'
           max_err_lower = MaxErrL'

return
