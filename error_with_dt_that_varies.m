clc, clear;
profile on
addpath('functions');
addpath('matrices');
addpath('lymph_minimal')

file_name = 'OutputData10p1.mat';
load(file_name);

T = 0.25;
Nt = 5;
inner_it = 2048;
eta = 10^(-0.5);
tol = 1e-3;
stagnation_tol = 0.3;

% Run the solver
[sigma_sol, ~, results] = CN_kaczmarz(file_name, T, Nt, 0.5, 'GRK', 1.4, tol, stagnation_tol, eta, 200, inner_it, 'N', 'Y');

dt = results.dt;
nts = size(results.num_quad_history, 2) - 1;
time_vector = 0:dt:(nts*dt);
time_vector_steps = dt:dt:(nts*dt); % For metrics evaluated only at the end of a step

% Initialize arrays for global relative errors over time
global_rel_error_gmres = zeros(1, length(time_vector));
global_rel_error_bs    = zeros(1, length(time_vector));

for k = 1:length(time_vector)
    % Extract the entire state vector for time step k
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

%%
% 1. Global Error Evolution
figure('Name', sprintf('Global Error (dt=%g)', dt), 'Position', [100, 200, 500, 400]);
hold on;
plot(time_vector, global_rel_error_gmres, '-o', 'LineWidth', 2, 'Color', '#A2142F', 'MarkerFaceColor', '#A2142F', 'DisplayName', 'GMRES Error');
plot(time_vector, global_rel_error_bs, '-s', 'LineWidth', 2, 'Color', '#77AC30', 'MarkerFaceColor', '#77AC30', 'DisplayName', 'Backslash Error');
grid on;
set(gca, 'YScale', 'log');
xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Global $L^2$ Relative Error', 'Interpreter', 'latex', 'FontSize', 12);
title('Error Evaluation (dt=%g)',dt,'FontSize', 13);
legend('Location', 'best');
% Print a quick summary to the command window
fprintf('\n--- Single Run Performance Summary ---\n');
fprintf('Max error: %4e\n', max(global_rel_error_gmres));
fprintf('Max error wrt backslash: %4e\n', max(norm(num_vec - bs_vec)));
fprintf('Total Inner Iterations: %d\n', sum(results.inner_iters));
fprintf('Total CPU Time: %.4f s\n', results.time_cpu_total);