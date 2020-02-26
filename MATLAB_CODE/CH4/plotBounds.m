 function [ Err ] = plotBounds( UBound, LBound, entries, diagentries, upbound, lowbound, Aparams, MaxErrU, t )
%PLOTBOUNDS plots the results of the obtained bounds.
%
%   function call:
%
% [ Err ] = plotBounds( UBound, LBound, entries, diagentries, upbound, lowbound )
%
%    input:
%
%         UBound: M x M matrix whose entries contain the upper bounds on the
%                 norms of the blocks of Z for a specific refinement step.
%         LBound: M x M diagonal matrix with the lower bounds on the norms
%                 of the blocks of Z in its entries for the specified
%                 refinement step.
%       lowbound: lower triangular matrix whose entries contain the upper
%                 bounds of the lower triangular off-diagonal blocks of Z.
%        upbound: upper triangular matrix whose entries contain the upper
%                 bounds of the upper triangular off-diagonal blocks of Z.
%        entries: M x M matrix whose entries contain the exact norms of the
%                 blocks of Z.
%    diagentries: 1 x M vector whose entries contain the exact norms of
%                 the diagonal blocks of Z.
%
%   output:
%
%            Err: R elative Error on the upper bound.
%
%            the function outputs some figures with the bounds (see below)
%
% Written by Carlos Echeverria. Last edited on May 4, 2018.


close all

 purple=[.49 .18 .55];
blue=[.301 .745 .933];

%% The following figure shows the upper and lower bounds on the diagonal elements:

figure(2),
           semilogy(upbound(1:end),'m:','LineWidth',4, 'Color', purple),hold on, grid on,
           semilogy(diagentries(1:end),'k-.','LineWidth',4),
           semilogy(lowbound(1:end),'b--','LineWidth',4, 'Color', blue),
           axis([1 9 1e-2 1e1])


 xlabel('block index i','FontSize',19),
 %ylabel('norm(Z_{i,i})','FontSize',19),
 %set(gca,'xtick',1:9),
 set(gca,'FontSize',20),
 leg = legend('upper bound (4.23)','norm(Z_{i,i})','lower bound (4.23)');
 set(leg,'FontSize',19,'Location','NorthEast' )
 title(sprintf('refinement step t=%d',t),'FontSize',20)

 if Aparams.Achoice==1, axis([1 length(upbound) 1e-2 1e2]); elseif Aparams.Achoice==2, axis([1 9 1e-3 1e-1]); elseif Aparams.Achoice==3, axis([1 9 1e-2 1e1]); elseif Aparams.Achoice==4, axis([1 9 1e-3 1e1]);end

    mkdir figures
    file_name = sprintf('/figures/9times9_Z%d_Bounds_t%d.eps',Aparams.Achoice, t);
    saveas(gcf,[pwd file_name],'epsc')
 %title('Bounds of block diagonal entries of A^{-1}','FontSize',16)
%  truesize(2,[300 300]);
 movegui('west')


%% The following figure shows the relative error of the upper bounds:

figure(3),
           Err = abs(UBound-entries)./UBound;
           surf(Err);

 xlabel('block index i','FontSize',19),
 ylabel('block index j','FontSize',19),
 zlabel('E_{ij}','FontSize',19),
 set(gca,'FontSize',20),
%leg = legend('Relative Error');
%set(leg,'FontSize',16,'Location','best' )
title(sprintf('refinement step t=%d',t),'FontSize',20)
    mkdir figures
    file_name = sprintf('/figures/9times9_Z%d_Error_t%d.eps',Aparams.Achoice, t);
    saveas(gcf,[pwd file_name],'epsc')
  movegui('east')

%% The following figure shows the maximum relative error of the upper bounds as a function of t

figure(4), clf

        semilogy(MaxErrU,':','LineWidth',5, 'Color', purple),

        grid on, axis([1 8 1e-15 1e0]),
        %xticks([1 2 3 4 5 6 7 8 9]),
        leg = legend('max_{i,j}(UB_{ij}-Z_{ij})/UB_{ij})');
        set(leg,'FontSize',16,'Location','West' )
        xlabel('refinement step t','FontSize',16),
        ylabel('Maximum relative error','FontSize',16)
        set(gca,'FontSize',20)
        title('Maximal Relative Eror on the Upper Bounds per refinement step','FontSize',16)
        movegui('south')




%% The following figure shows the upper and lower bounds on the diagonal
%  elements as well as the upper bound for the off-diagonal elements.
%  The image is a surface plot where the z-axis is the norm of the element
%  Z(i,j), while the x-axis is i-th row of the matrix Z and the y-axis is
%  the j-th column of the matrix Z:

% figure(5),
%            surf(UBound,'FaceColor',[1 0 0],'FaceAlpha',0.2), hold on,
%            surf(entries),
%            surf(LBound,'FaceColor',[0 1 0],'FaceAlpha',0.2); % alpha sets transparency
%
%  leg = legend('upper bound (2.23)','norm(Z_{ij})','lower bound (2.23)');
%  set(leg,'FontSize',16,'Location','best' )
%  xlabel('block index i','FontSize',16),
%  ylabel('block index j','FontSize',16),
%  zlabel('norm(Z_{i,i})','FontSize',16),
%  set(gca,'FontSize',20),
%  %title('Two-sided bounds on the block entries of A^{-1}','FontSize',16)
%  movegui('northwest')


%% The following figure shows the bounds on upper off diagonal elements.
%  To change the block row one needs to do it manually:

% figure(6),
%            semilogy(UBound(2,3:end),'r:','LineWidth',4),hold on, grid on,
%            semilogy(entries(2,3:end),'-.','LineWidth',4),
%
%  leg = legend('norm(Z_{2,j})', 'Bound (2.21)');
%  set(leg,'FontSize',16,'Location','best' )
%  xlabel('j-2','FontSize',16),
%  ylabel('norm(Z_{2,j})','FontSize',16),
%  set(gca,'FontSize',20),
%  %title('Row Decay of the norm of block entries of A^{-1}','FontSize',16)
%  movegui('southeast')
%

%% The following figure shows bounds on lower off diagonal elements
%  To change the block column one needs to do it manually:

% figure(7),
%           semilogy(UBound(3:end,2),'r:','LineWidth',4),hold on, grid on
%           semilogy(entries(3:end,2),'-.','LineWidth',4),
%
% leg = legend( 'Bound (2.22)', 'norm(Z_{i,2})');
% set(leg,'FontSize',16,'Location','southwest' )
% xlabel('i-2','FontSize',16),
% ylabel('norm(Z_{i,2})','FontSize',16),
% set(gca,'FontSize',20),
% %title('Column Decay of the norm of block entries of A^{-1}','FontSize',16)
% movegui('southwest')

end
