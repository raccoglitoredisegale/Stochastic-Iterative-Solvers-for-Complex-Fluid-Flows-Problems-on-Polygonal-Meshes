function [Uh, res_hist, adaptive_struct] = adaptive_F_AB_GMRES(A, b, u_0, choose_method, omega, tol, stagnation_tol, eta, m_max, initial_l_max, increment)
%
%  Adaptive Restart Flexible AB-GMRES method
%  (Dynamic inner-iteration tuning routine)
%   this method solves the linear system Au=b by executing successive restart cycles
%   of F_AB_GMRES. It monitors the convergence and stagnation profiles of the system. 
%   If progress stagnates above stagnation_tol, it scales the inner 
%   iteration limit (l_max) by the given increment factor for the next cycle.
%
%    INPUT:
%          A                : Matrix (Sparse format or not)
%          b                : vector
%          u_0              : initial solution
%          choose_method    : - K for (relaxed) Kaczmarz
%                             - GK for Greedy Kaczmarz
%                             - RK for Randomized Kaczmarz
%                             - GRK for Greedy Randomized Kaczmarz
%          omega            : relaxation parameter for Kaczmarz function
%          tol              : global tolerance for outer stopping criterion
%          stagnation_tol   : tolerance baseline to determine convergence stall
%          eta              : tolerance passed down to inner Kaczmarz loops
%          m_max            : maximum number of outer iterations per cycle
%          initial_l_max    : starting limit for inner iteration allocations
%          increment        : scaling factor to increase l_max upon stagnation
%
%    OUTPUT:
%          Uh               : final solution vector
%          res_hist         : accumulated residual history across all cycles
%          adaptive_struct  : structure containing execution statistics:
%                               - total_inner          : total inner updates
%                               - total_outer          : total outer iterations
%                               - stag_hist            : history of stagnation ratios
%                               - glob_conv_norm_vect  : relative residual history
%                               - max_inner_per_outer  : l_max tracking per cycle
%                               - global_res           : absolute residual tracking
%

Uh = u_0;
l_max = initial_l_max;
res_hist = [];

% Initialize  structure
adaptive_struct = struct();
adaptive_struct.total_inner = 0;
adaptive_struct.total_outer = 0;
adaptive_struct.stag_hist = [];
adaptive_struct.glob_conv_norm_vect = [];
adaptive_struct.max_inner_per_outer = [];
adaptive_struct.global_res = [];

r = b - A * Uh;
r_norm = norm(r);
norm_b = norm(b);
n_of_restart = 0;

% Restart sequence loop 
while adaptive_struct.total_inner < 12e6
    z_0 = zeros(size(Uh));
    
    % Call Flexible AB-GMRES with current parameter profile
    [z, ~, res, inner_iters, outer_iters] = F_AB_GMRES(A, r, z_0, choose_method, omega, tol, eta, m_max, l_max);
    
    % Update global performance statistics
    adaptive_struct.total_inner = adaptive_struct.total_inner + inner_iters;
    adaptive_struct.total_outer = adaptive_struct.total_outer + outer_iters;
    adaptive_struct.max_inner_per_outer = [adaptive_struct.max_inner_per_outer; l_max];
    res_hist = [res_hist; res];
    
    % Compute correction and evaluate new system state
    Uh = Uh + z;
    r = b - A * Uh;
    r_new_norm = norm(r);
    
    glob_conv_norm = r_new_norm / norm_b;
    stagnation_norm = r_new_norm / r_norm;
    
    adaptive_struct.stag_hist = [adaptive_struct.stag_hist; stagnation_norm];
    adaptive_struct.glob_conv_norm_vect = [adaptive_struct.glob_conv_norm_vect; glob_conv_norm];
    adaptive_struct.global_res = [adaptive_struct.global_res; r_new_norm];
    
    % Check global convergence target
    if glob_conv_norm < tol
        fprintf('Converged after %d restart cycles, with l_max=%d\n', n_of_restart, l_max);
        break;
    end
    
    % if stagnation: expand the inner preconditioning scope
    if stagnation_norm > stagnation_tol
        n_of_restart = n_of_restart + 1;
        l_max = round(increment * l_max);
    end
       
    r_norm = r_new_norm;
end
end