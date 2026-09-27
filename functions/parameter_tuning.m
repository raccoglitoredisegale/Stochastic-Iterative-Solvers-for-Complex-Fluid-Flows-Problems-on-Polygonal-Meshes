function [opt_l, opt_w] = parameter_tuning(A, b, eta, choose_method)
%
%  Parameter tuning for Kaczmarz-type methods
%  (algorithm 4.1)
%
%   this method finds the optimal number of inner row updates (opt_l) and relaxation 
%   parameter (opt_w) for a given Kaczmarz variant.
%
%    INPUT:
%          A                : Matrix (Sparse format or not)
%          b                : right hand side vector
%          eta              : tolerance for the inner stopping criterion
%          choose_method    : - K for (relaxed) Kaczmarz
%                             - GK for Greedy Kaczmarz
%                             - RK for Randomized Kaczmarz
%                             - GRK for Greedy Randomized Kaczmarz
%
%    OUTPUT:
%          opt_l            : optimal number of inner iterations (row updates)
%          opt_w            : optimal relaxation parameter omega
%

    [m, n] = size(A);
    norm_b = norm(b);
    row_norms_sq = sum(A.^2, 2);
    frob_norm_sq = sum(row_norms_sq);
    
    % Precompute residual update matrix
    C = A * A'; 
    
    % Randomized methods require multiple runs to find stable median values
    if strcmp(choose_method, 'RK') || strcmp(choose_method, 'GRK')
        num_runs = 10;
    else
        num_runs = 1;
    end
    
    l_max_list = zeros(num_runs, 1);
    w_opt_list = zeros(num_runs, 1);
    w_candidates = 0.05:0.05:1.5;
    
    for run = 1:num_runs
        % Phase 1: Estimate l_max using baseline omega = 1.0
        l = 0;
        res = b; 
        w_temp = 1.0;
        
        while norm(res) > eta * norm_b
            l = l + 1;
            
            switch choose_method
                case 'K'
                    i = mod(l - 1, m) + 1;
                case 'GK'
                    [~, i] = max(abs(res));    
                case 'RK'
                    prob_dist = row_norms_sq / frob_norm_sq;
                    cdf = cumsum(prob_dist); 
                    i = find(cdf >= rand, 1);
                case 'GRK' 
                    res_norm_sq = norm(res)^2;
                    epsilon_p = 0.5 * ((1/res_norm_sq) * ...
                        max(abs(res).^2 ./ row_norms_sq) + (1/frob_norm_sq));
                    U_p_indices = find(abs(res).^2 >= epsilon_p * res_norm_sq * row_norms_sq);
                    
                    if isempty(U_p_indices)
                        [~, i] = max(abs(res)); 
                    else
                        weights = abs(res(U_p_indices)).^2;
                        prob_dist = weights / sum(weights);
                        cdf = cumsum(prob_dist);
                        selected_idx = find(cdf >= rand, 1);
                        i = U_p_indices(selected_idx);
                    end
            end
            
            % Update residual vector without full matrix re-multiplication
            res = res - w_temp * (res(i)/row_norms_sq(i)) * C(:, i);
            
            if l > 50 * m 
                break;
            end
        end
        l_max_list(run) = l;
        
        % Phase 2: Search for optimal omega using the fixed l_max profile
        best_res = inf;
        best_w = inf;
        
        for w = w_candidates
            res = b; 
            
            for p = 1:l
                switch choose_method
                    case 'K'
                        i = mod(p - 1, m) + 1;
                    case 'GK'
                        [~, i] = max(abs(res));    
                    case 'RK'
                        prob_dist = row_norms_sq / frob_norm_sq;
                        cdf = cumsum(prob_dist); 
                        i = find(cdf >= rand, 1);
                    case 'GRK' 
                        res_norm_sq = norm(res)^2;
                        epsilon_p = 0.5 * ((1/res_norm_sq) * ...
                            max(abs(res).^2 ./ row_norms_sq) + (1/frob_norm_sq));
                        U_p_indices = find(abs(res).^2 >= epsilon_p * res_norm_sq * row_norms_sq);
                        if isempty(U_p_indices)
                            [~, i] = max(abs(res));
                        else
                            weights = abs(res(U_p_indices)).^2;
                            prob_dist = weights / sum(weights);
                            cdf = cumsum(prob_dist);
                            selected_idx = find(cdf >= rand, 1);
                            i = U_p_indices(selected_idx);
                        end
                end
                res = res - w * (res(i)/row_norms_sq(i)) * C(:, i);
            end
            
            if norm(res) < best_res
                best_res = norm(res);
                best_w = w;
            end
        end
        w_opt_list(run) = best_w;
    end
    
    % Compute median results across sample runs
    opt_l = round(median(l_max_list));
    opt_w = median(w_opt_list);
    
    % Standard Kaczmarz sweeps full cycles; round up to nearest complete pass
    if strcmp(choose_method, 'K')
        opt_l = ceil(opt_l / m) * m; 
    end
    
    fprintf('Tuning [%s]: opt. inner it = %d, opt. w = %.1f\n', choose_method, opt_l, opt_w);
end