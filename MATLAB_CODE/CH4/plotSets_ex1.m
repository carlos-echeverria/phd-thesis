function [  ] = plotSets_ex1( X, Y, FV, ELN, d)
%PLOTSETS Plots the inlcusion sets for the specific example given in
% section 3 of the paper:
%
%  [ Echeverria, Liesen, Nabben - Block diagonal dominance of matrices
%    revisited: bounds for the norms of inverses and eigenvalue
%    inclusion sets - Linear Algebra and its applications - 2018 ]
%
%   function call:
%
%      [  ] = plotSets(X, Y, FVV, ELNN1, ELNN2, d, e)
%
%    input:
%
%          X: X-coordinates of a rectangular grid (X, Y)
%          Y: Y-coordinates of a rectangular grid (X, Y)
%        FVV: the function of the complex plane evaluated at the points z=x+iy
%       ELNN: the function of the complex plane evaluated at the points z=x+iy
%          d: vector with the eigenvalues of A
%          e: vector of zeros of size(d).
%
%   output:
%
%          figures with the inclusion sets.
%
% Written by Carlos Echeverria. Last edited on October 25, 2019.


figure(8),clf,
purple = [.49 .18 .55];
blue = [.301 .745 .933];

evals=d;
p1=plot(real(d), imag(d), 'k.', 'MarkerSize', 28); hold on,
linS = {'--',':','-.','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};


    [C1,h1] = contour(X, Y,  FV{1},  [1 1], 'LineStyle', linS{1}, 'LineWidth', 3, 'Color', purple);
    [C1,h1] = contour(X, Y,  FV{2},  [1 1], 'LineStyle', linS{2}, 'LineWidth', 3, 'Color', purple);

    %clabel(C,h,'FontSize',10);

axis equal, grid on,
leg = legend('eigenvalues', 'G_{1}^{FV}', 'G_{2}^{FV}');
set(leg, 'FontSize', 17, 'Location', 'NorthOutside', 'Orientation', 'horizontal');
legend boxoff;

xlabel('Re(z)','FontSize',17),
ylabel('Im(z)','FontSize',17),
%title('Eigenvalue Inclusion Regions','FontSize',17)
set(gca,'FontSize',18),
% movegui('northeast')




figure(9),clf,
purple = [.49 .18 .55];
blue = [.301 .745 .933];

p1=plot(real(d), imag(d), 'k.', 'MarkerSize', 28); hold on,
linS = {'--',':','-.','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};
% linS = fliplr(linS);

%     txt = sprintf('G_{%d}^{ELN}', ii);
    [C2,h2] = contour(X, Y,  ELN{1},  [1 1], 'LineStyle', linS{1}, 'LineWidth', 3,'Color', blue);
    [C2,h2] = contour(X, Y,  ELN{2},  [1 1], 'LineStyle', linS{2}, 'LineWidth', 3,'Color', blue);

    %clabel(C,h,'FontSize',10);



axis equal, grid on,
leg = legend('eigenvalues','G_{1}^{new}', 'G_{2}^{new}');
set(leg, 'FontSize', 17, 'Location', 'NorthOutside', 'Orientation', 'horizontal');
legend boxoff;

xlabel('Re(z)','FontSize',17),
ylabel('Im(z)','FontSize',17),
%title('Eigenvalue Inclusion Regions','FontSize',17)
set(gca,'FontSize',18),
% movegui('northeast')





% figure(9),clf,
%         surfc(X,Y,ELN{2});
%
%         xlabel('Re(z)','FontSize',16),
%         ylabel('Im(z)','FontSize',16),zlabel('f(z)','FontSize',16)
%         title('f(z)=||A_{ij}|| \cdot ||(A_{ii}-zI)^{-1}||','FontSize',16)
        % movegui('northeast')

% figure(10),clf,
%       contour(X,Y,FVV,[1 1],'r'),
%
%        axis equal;
%        leg = legend('||A_{ij}|| \cdot ||(A_{ii}-zI)^{-1}||=1');
%        set(leg,'FontSize',12,'Location','northwest' ),
%        xlabel('Re(z)','FontSize',16),
%        ylabel('Im(z)','FontSize',16),
%        title('Eigenvalue Inclusion Regions','FontSize',16),

% figure(11),clf,
%       surfc(X,Y,ELNN);
%
%        xlabel('Re(z)','FontSize',16),
%        ylabel('Im(z)','FontSize',16),
%        zlabel('f(z)','FontSize',16),
%        title('f(z)=||(A_{ii}-zI)^{-1}\cdotA_{ij}||','FontSize',16),
%        movegui('southwest')
%

% figure(12),clf,
%       contour(X,Y,ELNN,[1 1],'b-.'),
%
%        axis equal,
%        leg = legend('||A_{ij}\cdot(A_{ii}-zI)^{-1}||=1','||(A_{ii}-zI)^{-1}\cdotA_{ij}||=1');
%        set(leg,'FontSize',12,'Location','northwest' )
%        xlabel('Re(z)','FontSize',16),
%        ylabel('Im(z)','FontSize',16)
%        title('Eigenvalue Inclusion Regions','FontSize',16)

end
