clc; clear; 
profile on
addpath('functions');
addpath('matrices');
addpath('lymph_minimal')


% --- Setup Parameters ---
file_name = 'OutputData10p1.mat';
load(file_name);

T = 0.05;
Nt = 1;
eta = 10^(-0.5);
tol = 1e-3;
stagnation_tol = 0.3;


% Array of initial inner iterations to test
inner_it_array = [64, 128, 192, 256, 384, 512, 768, 1024, 2048, 4096, 8192, 15000];
num_tests = length(inner_it_array);
    
% Arrays to store the metrics for the graphs
time_taken = zeros(1, num_tests);
total_inner_iterations = zeros(1, num_tests);
max_error_gmres_array = zeros(1, num_tests);
max_error_bs_array = zeros(1, num_tests);

fprintf('\n=== Starting Adaptive GMRES Analysis ===\n');

for i = 1:num_tests
    inner_it = inner_it_array(i);
    fprintf('\n--- Testing with initial inner iterations = %d ---\n', inner_it);
    
    % Run the solver
    [sigma_sol, ~, results] = CN_kaczmarz(file_name, T, Nt, 0.5, 'GK', 1.4, tol, stagnation_tol, eta, 200, inner_it, 'N', 'Y');
    
    % Reconstruct error evaluation for this specific run
    nts = size(results.num_quad_history, 2) - 1;
    global_rel_error_gmres = zeros(1, nts + 1);
    global_rel_error_bs    = zeros(1, nts + 1);
    
    for k = 1:(nts + 1)
        num_vec = results.num_quad_history(:, k);
        bs_vec  = results.bs_quad_history(:, k);
        ex_vec  = results.exact_quad_history(:, k);
        
        n_ex = norm(ex_vec);
        
        if n_ex < 1e-14
            global_rel_error_gmres(k) = norm(num_vec - ex_vec);
            global_rel_error_bs(k)    = norm(bs_vec - ex_vec);
        else
            global_rel_error_gmres(k) = norm(num_vec - ex_vec) / n_ex;
            global_rel_error_bs(k)    = norm(bs_vec - ex_vec) / n_ex;
        end
    end
    
    % Extract the max errors for the current run
    current_max_gmres = max(global_rel_error_gmres);
    current_max_bs    = max(global_rel_error_bs);
    
    % Store data for plotting
    time_taken(i) = results.time_cpu_total;
    total_inner_iterations(i) = sum(results.inner_iters);
    max_error_gmres_array(i) = current_max_gmres;
    max_error_bs_array(i) = current_max_bs;
    
    % Console output
    fprintf('Max error: %e\n', current_max_gmres);
    fprintf('Max error backslash: %e\n', current_max_bs);
end

fprintf('\n=== Analysis Complete, Generating Plots ===\n');

% --- Plotting the 3 Separate Graphs ---

% Graph 1: Time vs Max Inner Iterations
figure('Name', 'CPU Time Analysis', 'Position', [100, 200, 500, 400]);
plot(inner_it_array, time_taken, '-o', 'LineWidth', 5, 'Color', '#A2142F', 'MarkerFaceColor', '#A2142F');
set(gca, 'FontSize', 20);
set(gca, 'XScale', 'log'); 
grid on;
xlabel('Initial Inner Iterations ($l_{max}$)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Total CPU Time (s)', 'Interpreter', 'latex', 'FontSize', 20);
title('CPU Time vs Initial $l_{max}$', 'Interpreter', 'latex', 'FontSize', 20);

% Graph 2: Total Inner Iterations vs Max Inner Iterations
figure('Name', 'Inner Iterations Analysis', 'Position', [620, 200, 500, 400]);
plot(inner_it_array, total_inner_iterations, '-s', 'LineWidth', 5, 'Color', '#0072BD', 'MarkerFaceColor', '#0072BD');
set(gca, 'XScale', 'log');
grid on;
xlabel('Initial Inner Iterations ($l_{max}$)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Total Inner Iterations Performed', 'Interpreter', 'latex', 'FontSize', 20);
title('Total Inner Iters vs Initial $l_{max}$', 'Interpreter', 'latex', 'FontSize', 20);

% Graph 3: Max Error vs Max Inner Iterations
figure('Name', 'Solver Error Analysis', 'Position', [1140, 200, 500, 400]);
loglog(inner_it_array, max_error_gmres_array, '-^', 'LineWidth', 5, 'Color', '#D95319', 'MarkerFaceColor', '#D95319', 'DisplayName', 'GMRES Max Error');
hold on;
loglog(inner_it_array, max_error_bs_array, '--k', 'LineWidth', 5, 'DisplayName', 'Backslash Max Error');
grid on;
xlabel('Initial Inner Iterations ($l_{max}$)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Max Global Relative Error', 'Interpreter', 'latex', 'FontSize', 20);
title('Solver Error vs Initial $l_{max}$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Location', 'best', 'FontSize', 20);