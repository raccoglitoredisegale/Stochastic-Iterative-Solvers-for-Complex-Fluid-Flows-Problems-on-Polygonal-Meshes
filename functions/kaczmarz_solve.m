function [z_k, res_k, l, z_history, res_hist] = kaczmark_solve (A, v_k, l_max, eta, choose_method, omega)
%
%  Kaczmarz-type inner-iteration algorithms
%   this function computes a preconditioned direction vector z_k by solving 
%   the inner linear subproblem Az = v_k. It supports standard cyclic, greedy, 
%   randomized, or greedy randomized row projection methods.
%
%    INPUT:
%          A                : Matrix (Sparse format or not)
%          v_k              : right-hand side vector (current GMRES basis direction)
%          l_max            : maximum number of inner iterations (row updates)
%          eta              : inner loop relative tolerance limit
%          choose_method    : - K for (relaxed) Kaczmarz
%                             - GK for Greedy Kaczmarz
%                             - RK for Randomized Kaczmarz
%                             - GRK for Greedy Randomized Kaczmarz
%          omega            : relaxation parameter
%
%    OUTPUT:
%          z_k              : final preconditioned direction vector
%          res_k            : absolute final residual norm
%          l                : total number of inner iterations performed
%          z_history        : matrix tracking the solution vector at each projection step
%          res_hist         : array tracking relative residual norm variations
%

[m,n] = size(A);
z_k = zeros(n,1);                   % initial approximate solution
l = 0;                              % number of iterations
z_history = zeros(n, l_max); 
row_norms_sq = sum(A.^2, 2);        % squared norms of rows ||alpha_i||^2
frob_norm_sq = sum(row_norms_sq);   % Squared Frobenius norm ||A||^2
res = v_k - A*z_k;                  % vector of residuals (one per row of A)
res_k = norm(res,2);

% Preallocate residual history
res_hist = zeros(l_max+1, 1);
res_hist(1) = res_k;

C = A*A';                           % for calculating the residual in the end
% the algorithm keeps running until either:
%   - we reach the maximum number of iteration
%   - the residual error is eta*||v_k|| 
while (l < l_max && norm(res)/norm(v_k) >= eta)
    l = l+1;
    % fprintf('   inner iter: %d \n',l);
    switch choose_method
        case 'K' % (relaxed) Kaczmarz: cycle through every row
            i = mod(l - 1, m) + 1;
        case 'GK'
            % Greedy Kaczmarz: pick the row that maximizes the residual
            [~, i] = max(abs(res));
        case 'RK'
            % Randomized Kaczmarz: pick rows based on a probability
            % equal to ||row||^2 / ||A||^2
            prob_dist = row_norms_sq / frob_norm_sq;
            cdf = cumsum(prob_dist); % cumulative distribution
            i = find(cdf >= rand, 1);
        case 'GRK' % Greedy Randomized Kaczmarz
            res_norm_sq = norm(res)^2;
            % threshold:
            epsilon_p = 0.5 * ((1/res_norm_sq) * ...
                max(abs(res).^2 ./ row_norms_sq) + (1/frob_norm_sq));
            % filtering out the rows where the residual is not worth
            % considering:
            U_p_indices = find(abs(res).^2 >= epsilon_p * res_norm_sq * row_norms_sq);
            % for the probabilities: only considering rows in U_p_indices
            weights = abs(res(U_p_indices)).^2;
            prob_dist = weights / sum(weights);
            %from now on analogous as RK case
            cdf = cumsum(prob_dist);
            selected_idx = find(cdf >= rand, 1);
            i = U_p_indices(selected_idx);
        otherwise
            error('Unknown method selected.');
    end
    % common to all methods:
    z_k = z_k + omega * (res(i)/row_norms_sq(i)) * (A(i,:))';
    
    z_history(:, l) = z_k;
    
      res = res - omega * (res(i)/row_norms_sq(i)) * C(:, i);
    
    % Store the relative Euclidean L2 norm for the GMRES termination criteria tracking
    res_hist(l+1) = norm(res, 2)/norm(v_k);
end

res_k = norm(res, 2); 

% if the solver stops early we remove the preallocated space
z_history = z_history(:, 1:l);
res_hist = res_hist(1:l+1);
end