% Written by Carlos Echeverria on December 9, 2019.

prob = 1;       %  1: Conv-Diff
   N = 40;      %  number of intervals of the grid in x-dir (must be even)
   M = 40;      %  number of intervals of the grid in y-dir (must be even)
 problem_size = (N-1)*(M-1)

  epsilons = [1e-8, 1e-6, 1e-4]; % desired perturbation parameters
local_tols = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-13];
% local_tols = [1e-2];


for i = 1:length(epsilons)
    epsi = epsilons(i);

    experiments_inexact_mSm_upwind(prob, N, M, epsi, local_tols);
    % max_err_vec(j)=max_err;
    file_name = sprintf('mSm_inexact_eps_%5.0e_N_%d_M_%d.eps', epsi, N, M);
    saveas(gcf, file_name,'epsc')

end
