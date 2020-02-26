function [ ] = block_Schwarz_ploterr_inexact(EQparams, SOLUTIONparams)
%BLOCK_SCHWARZ_PLOTERR plots the convergence history of errors generated
% by the Multiplicative Schwarz method as well as the theoretical error
% bounds for both orderings of the iteration matrix for each step of the
% iteration.
%
%
%   function call:
%
%       [ ] = block_Schwarz_ploterr(EQparams, SOLUTIONparams)
%
%   input:
%
%        EQparams: data structure the following equation parameters:
%        problem: problem identifyer (1 for conv-diff, 2 for Laplace)
%
%  SOLUTIONparams: data structure with the following solution parameters:
%         err_s1: vector with the computed error of the mSm at each step
%                 for T=Q2Q1.
%         err_s2: vector with the computed error of the mSm at each step
%                 for T=Q1Q2.
%   err_bound_s1: vector with the algebraic error bound of the mSm at each
%                 step for T=Q1Q2.
%   err_bound_s2: vector with the algebraic error bound of the mSm at each
%                 step for T=Q2Q1.
%
%   output:
%
%    figure plot showing the error and its corresponding bound at each
%    iteration step of the multiplicative Schwarz method.
%
%
% Written by Carlos Echeverria on September 23, 2016.
% Last Edited by C. E. on September 24, 2019.

  prob = EQparams.problem;
err_s1 = SOLUTIONparams.err_s1;
err_s2 = SOLUTIONparams.err_s2;

figure(101), clf;
semilogy(0:length(err_s2)-1, err_s2/err_s2(1), 'k-.', 'LineWidth', 2); hold on;
semilogy(0:length(err_s1)-1, err_s1/err_s1(1), 'k-' , 'LineWidth', 2);

   leg = legend('inexact', 'T = Q_2Q_1', 'Location', 'Best');
   set(leg,'FontSize',17);
   set(gca,'FontSize',18);
xlabel('k','FontSize',17);
 title('Upwind','FontSize',17);
% axis([0 20 1e-16 1e0]);
set(gca,'XTick',[0:1:10,12:2:16]')

if prob==1
    ylabel('Error norms and bounds (||\cdot||_{\infty})','FontSize',17);
elseif prob == 2
    ylabel('Error norms and bounds (||\cdot||_{2})','FontSize',17);
end

%titletext = ['\epsilon = ',num2str(epsi)];
%title(titletext,'FontSize',16);


end
