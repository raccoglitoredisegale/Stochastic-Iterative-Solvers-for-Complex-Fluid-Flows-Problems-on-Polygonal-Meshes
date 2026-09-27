clc, clear;

addpath('functions');
addpath('matrices');
addpath('lymph_minimal')

% Load data
data_file = 'OutputData10p1.mat';
load(data_file);
T = 0.05;
Nt = 1;
inner_it_vect = [128, 256, 512, 768, 1024, 1536, 2048, 3072, 4096, 8192, 16384, 32768, 60000, 90000];


% Assign the figure to a variable 'fig' so it can be referenced for saving
fig = figure('Name', 'Preconditioner Efficiency', 'Color', 'w', 'Position', [100, 100, 1920, 1080]);
hold on;
grid on;
colororder(turbo(length(inner_it_vect)));

for i = 1:length(inner_it_vect)
    inner_it = inner_it_vect(i);
    fprintf('\n--- Running Test %d/%d (inner_it = %d) ---\n', i, length(inner_it_vect), inner_it);
    
    [~, ~, results] = CN_kaczmarz(data_file, T, Nt, 0.5, 'GRK', 1.4, 1e-3, 0, 10^(-0.5), 200, inner_it, 'N', 'N');
    
    res_history = results.res_history{1};
    rel_res = res_history / res_history(1);
    iterations = 0:(length(rel_res)-1);
    
    plot(iterations, rel_res, 'LineWidth', 5, 'MarkerSize', 5, ...
        'DisplayName', sprintf('l_{max} = %d', inner_it));
end

set(gca, 'YScale', 'log', 'FontSize', 20);
xlabel('Outer GMRES Iterations', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Relative Residual $$||b - Ax_k||_2 / ||r_0||_2$$', 'Interpreter', 'latex', 'FontSize', 20);
title(sprintf('F-AB-GMRES Convergence (GRK Preconditionerz)'), 'FontSize', 20);
legend('Location', 'northeast', 'FontSize', 20);
hold off;