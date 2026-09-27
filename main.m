clc, clear; 
addpath('functions');

tol = 1e-6;
outer_it_max = 1000;
eta = 10^(-2);
matrix_name = 'RANDL1';

A = mmread(['./matrices/500x5000/', matrix_name, '.mtx']);
[m, n] = size(A);

x_true = rand(n, 1);
b = A * x_true;
u0 = rand(n, 1);

cond = condition_number(A);

fprintf("Eta: %.1e\n", eta);

methods = {'NESOR','GRK','GK'};
num_methods = length(methods);

results_res = cell(num_methods, 1);
stats = struct();

tic;
C = A * A';
time_C = toc;

for j = 1:num_methods
    meth = methods{j};
    tic;
    if strcmp(meth, 'NESOR')
        [opt_l, opt_w] = parameter_tuning(A, b, eta, 'K');
        tune_time = toc;
        inner_it_max = max(1, ceil(opt_l / m));

        tic;
        [u_out, res, total_inner,avg_inner] = AB_GMRES(A, b, u0, opt_w, tol, outer_it_max, inner_it_max);
        solve_time = toc;
    else
        %[opt_l, opt_w] = parameter_tuning(A, b, eta, meth)
        tune_time = toc;
        opt_w = 1.0;
        tic;
        [u_out, ~, res, total_inner,~,avg_inner] = F_AB_GMRES(A, b, u0, meth, 1.0, tol, eta, outer_it_max, 10);
        solve_time = toc;
    end

    stats.(meth).name = meth;
    stats.(meth).outer = length(res) - 1;
    stats.(meth).inner = total_inner;
    stats.(meth).w = opt_w;
    stats.(meth).tune = tune_time;
    stats.(meth).solve = solve_time;
    stats.(meth).total = tune_time + solve_time + time_C;
    stats.(meth).rel_res = norm(b - A*u_out) / norm(b);
    stats.(meth).rel_err = norm(u_out - x_true) / norm(x_true);
    stats.(meth).avg = avg_inner;

    results_res{j} = res;
end

figure('Color', 'w');
hold on;
colors = lines(num_methods);
line_styles = {'-', '-', ':', '-.', '--'};

for j = 1:num_methods
    res = results_res{j};
    num_iters = length(res);
    semilogy(0:num_iters-1, res/res(1), ...
        'Color', colors(j,:), ...
        'LineStyle', line_styles{j}, ...
        'LineWidth', 2.0, ...
        'DisplayName', methods{j});
end

grid on;
set(gca,'YScale','log');
set(gca, 'GridAlpha', 0.25, 'MinorGridAlpha', 0.15, 'FontSize', 12);
xlabel('Number of Outer iterations $k$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Relative residual norm $\|r_k\| / \|r_0\|$', 'Interpreter', 'latex', 'FontSize', 13);
title('Convergence of preconditioned AB-GMRES methods', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'Interpreter', 'latex');
hold off;

fprintf('\n\n');
fprintf('%-15s %-15s | %s\n', 'Outer', 'Inner', matrix_name);
fprintf('%-15s %-15s | %.2f\n', 'iteration', 'iteration', cond);
fprintf('------------------------------------------------------\n');
for j = 1:num_methods
    meth = methods{j};
    safe_meth = strrep(meth, '-', '');
    s = stats.(safe_meth);
    
    % Format the first column string
    if strcmp(meth, 'NE-SOR')
        col1 = 'AB-GMRES';
    elseif j == 2
        col1 = 'F-AB-GMRES'; % Print algorithm name only on the first Kaczmarz row
    else
        col1 = '';
    end
    
    fprintf('%-15s %-15s | %-5d (%d, %d, %.1f)\n', col1, meth, s.outer, s.inner, s.avg, s.w);
    fprintf('%-15s %-15s | %-5.2f      (%.2f)\n', '', '', s.total, s.tune);
    fprintf('%-15s %-15s | %15.2e      %15.2e\n \n', '', '', s.rel_res, s.rel_err);
    
    if j < num_methods
        fprintf('------------------------------------------------------\n');
    end
end
fprintf('------------------------------------------------------\n');


%%
clc;clear;
addpath('functions');
addpath('matrices');


file_name = 'OutputData10p1.mat';
load(file_name);

T = 0.25;
Nt = 5;
inner_it = 1024;
eta = 1e-1;
tol = 1e-6;
stagnation_tol = 0.30;
                        % CN_kaczmarz(data_file, T, Nt, theta, choose_method, omega, tol, eta, outer_it_max, inner_it_max, PT, adaptive?)
[sigma_sol, ~, results] = CN_kaczmarz(file_name, T, Nt, 0.5, 'GK', 1.1, tol, stagnation_tol, eta, 200, inner_it, 'Y', 'N');

% to check results for any point in the mesh:
CN_plot(results, 10); 


%% Kaczmarz analysis dt

addpath('functions');
addpath('matrices');
data_file = 'OutputData10p1.mat'; 
T = 0.25;

clear L2_error MAX_error avg_outer avg_inner total_inner_cost total_outer_cost c_vect;

% choose dt:
Nt_vec = [5,10,15];
dt_vec = T./Nt_vec;

% F-AB-GMRES parameters
choose_method = 'GK'; 
omega = 1.0; 
tol = 1e-5; 
stagnation_tol = 0.3; 
eta = 10^(-0.5); 
outer_it_max = 200; 
inner_it_max = 1024;
Parameter_Tuning = 'Y';
adaptive = 'N'; % Default execution strategy choice

fprintf('Starting dt convergence analysis...\n');

for i = 1:length(Nt_vec)
    Nt_current = Nt_vec(i);
    fprintf('\n--- Running Test %d/%d (dt = %g) ---\n', i, length(Nt_vec), dt_vec(i));
    
    % Run the core solver matching the updated 13-input parameter signature
    [~, ~, res] = CN_kaczmarz(data_file, T, Nt_current, 0.5, choose_method, ...
                        omega, tol, stagnation_tol, eta, outer_it_max, inner_it_max, ...
                        Parameter_Tuning, adaptive);
    
    err_history = res.global_rel_error;
    err_history(isnan(err_history) | isinf(err_history)) = 0; 
    
    % --- DATA COLLECTION ---
    L2_error(i)  = rms(err_history); % Global performance metric
    MAX_error(i) = max(err_history); % Worst-case spike metric
    
    % iterations metrics
    c_vect(i)           = res.c;
    avg_outer(i)        = mean(res.outer_iters);
    avg_inner(i)        = mean(res.inner_avg);
    total_outer_cost(i) = sum(res.outer_iters);
    total_inner_cost(i) = sum(res.inner_iters); % Total inner iterations across the whole run
end

% making all vectors vertical
dt_plot       = dt_vec(:);
L2_plot       = L2_error(:);
MAX_plot      = MAX_error(:);
inner_total   = total_inner_cost(:);
outer_total   = total_outer_cost(:);
inner_to_outer_ratio = inner_total ./ outer_total;
outer_average = avg_outer(:);
inner_average = avg_inner(:);
c_vect        = c_vect(:);

%% --- PLOTS --- 

% 1. CONVERGENCE OF ERROR (LOGLOG)
figure('Name', 'Error Convergence Analysis', 'Position', [100, 100, 550, 400]);
loglog(dt_plot, L2_plot, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', ' $L^2$ Error');
hold on;
loglog(dt_plot, MAX_plot, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Max Error ($L^\infty$)');
loglog(dt_plot, (dt_plot.^2) * (L2_plot(1) / dt_plot(1)^2), 'k:', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{O}(\Delta t^2)$ Slope');
grid on;
xlabel('$\Delta t$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Relative Error', 'FontSize', 12);
title('F-AB-GMRES Time Convergence', 'Interpreter', 'latex', 'FontSize', 13);
legend('Location', 'best', 'Interpreter', 'latex');

% Grab standard system color matrix for the rest of the figures
colors = get(gca, 'ColorOrder'); 

% 2. CONDITIONING NUMBER AT EACH dt (NORMAL PLOT - Color Index 1)
figure('Name', 'System Matrix Conditioning', 'Position', [150, 150, 550, 400]);
plot(dt_plot, c_vect, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', colors(1,:), 'MarkerFaceColor', colors(1,:));
grid on;
set(gca, 'XDir', 'reverse');
xlabel('$\Delta t$ (Refining $\rightarrow$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Conditioning Number', 'FontSize', 12);
title('Matrix Condition Number vs $\Delta t$', 'Interpreter', 'latex', 'FontSize', 13);

% 3. TOTAL OUTER ITERATIONS AT EACH dt (NORMAL PLOT - Color Index 2)
figure('Name', 'Total Outer Iterations', 'Position', [200, 200, 550, 400]);
plot(dt_plot, outer_total, '-s', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', colors(2,:), 'MarkerFaceColor', colors(2,:));
grid on;
set(gca, 'XDir', 'reverse');
xlabel('$\Delta t$ (Refining $\rightarrow$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Total Outer Iterations performed', 'FontSize', 12);
title('Total Outer GMRES Iterations vs $\Delta t$', 'Interpreter', 'latex', 'FontSize', 13);

% 4. TOTAL INNER/TOTAL OUTER AT EACH dt (NORMAL PLOT - Color Index 3)
figure('Name', 'Inner over Outer Ratio', 'Position', [250, 250, 550, 400]);
plot(dt_plot, inner_to_outer_ratio, '-d', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', colors(3,:), 'MarkerFaceColor', colors(3,:));
grid on;
set(gca, 'XDir', 'reverse');
xlabel('$\Delta t$ (Refining $\rightarrow$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Total Inner / Total Outer Iterations', 'FontSize', 12);
title('Average Work Ratio per Outer Dimension step', 'Interpreter', 'latex', 'FontSize', 13);

fprintf('\n======================================= SUMMARY STATISTICS TABLE =======================================\n');
SummaryTable = table(dt_plot, c_vect, outer_total, outer_average, inner_total, inner_average, inner_to_outer_ratio, L2_plot, ...
    'VariableNames', {'dt', 'Cond_Number', 'Total_Outer_Its', 'Avg_Outer_Per_Step', 'Total_Inner_Its', 'Avg_Inner_Per_Step', 'Inner_Outer_Ratio', 'L2_Error'});
disp(SummaryTable);