function CN_plot(results, target_pt)
% CN_PLOT Generates analysis plots for Crank-Nicolson F-AB-GMRES
%
%   INPUT:
%       results   : Output struct from CN_kaczmarz
%       target_pt : (Optional) Quadrature point index to track. Default is 10.

    if nargin < 2
        target_pt = 10;
    end

    T=results.T;
    dt=results.dt;
    num_quad_points = size(results.num_quad_history, 1) / 4;
    nts = size(results.num_quad_history, 2) - 1;
    time_vector = 0:dt:(nts*dt);

    % check if target is within bounds
    if target_pt > num_quad_points || target_pt < 1
        fprintf('Warning: target_pt %d is out of bounds. Defaulting to middle point.\n', target_pt);
        target_pt = round(num_quad_points / 2); 
    end

    % row indices for target point
    idx_s11 = target_pt;
    idx_s12 = target_pt + num_quad_points;
    idx_s21 = target_pt + 2 * num_quad_points;
    idx_s22 = target_pt + 3 * num_quad_points;
    
    % GMRES tracking lines
    num_s11 = results.num_quad_history(idx_s11, :);
    num_s12 = results.num_quad_history(idx_s12, :);
    num_s21 = results.num_quad_history(idx_s21, :);
    num_s22 = results.num_quad_history(idx_s22, :);
    
    % Backslash tracking lines
    bs_s11 = results.bs_quad_history(idx_s11, :);
    bs_s12 = results.bs_quad_history(idx_s12, :);
    bs_s21 = results.bs_quad_history(idx_s21, :);
    bs_s22 = results.bs_quad_history(idx_s22, :);
    
    % Exact slution otracking lines
    ex_s11 = results.exact_quad_history(idx_s11, :);
    ex_s12 = results.exact_quad_history(idx_s12, :);
    ex_s21 = results.exact_quad_history(idx_s21, :);
    ex_s22 = results.exact_quad_history(idx_s22, :);
    
    % relative errors over time
    pt_rel_error_gmres = zeros(1, length(time_vector));
    pt_rel_error_bs    = zeros(1, length(time_vector));
    
    for k = 1:length(time_vector)
        num_vec = [num_s11(k); num_s12(k); num_s21(k); num_s22(k)];
        bs_vec  = [bs_s11(k);  bs_s12(k);  bs_s21(k);  bs_s22(k)];
        ex_vec  = [ex_s11(k);  ex_s12(k);  ex_s21(k);  ex_s22(k)];
        
        n_ex = norm(ex_vec);
        pt_rel_error_gmres(k) = norm(num_vec - ex_vec) / n_ex;
        pt_rel_error_bs(k)    = norm(bs_vec - ex_vec) / n_ex;
 
    end
    
    % --- Plotting ---
    figure('Name', sprintf('Point %d (dt=%f)', target_pt, dt), 'Position', [100, 100, 1200, 500]);
    
    % Subplot 1: Evolution of all 4 components
    subplot(1, 2, 1);
    hold on;
    
    % GMRES (Solid Lines)
    plot(time_vector, num_s11, '-', 'LineWidth', 2, 'Color', '#0072BD', 'DisplayName', 'GMRES \sigma_{11}');
    plot(time_vector, num_s12, '-', 'LineWidth', 2, 'Color', '#D95319', 'DisplayName', 'GMRES \sigma_{12}');
    plot(time_vector, num_s21, '-', 'LineWidth', 2, 'Color', '#EDB120', 'DisplayName', 'GMRES \sigma_{21}');
    plot(time_vector, num_s22, '-', 'LineWidth', 2, 'Color', '#7E2F8E', 'DisplayName', 'GMRES \sigma_{22}');
    
    % Backslash (Dotted Lines)
    plot(time_vector, bs_s11, ':', 'LineWidth', 2.5, 'Color', '#0072BD', 'HandleVisibility', 'off');
    plot(time_vector, bs_s12, ':', 'LineWidth', 2.5, 'Color', '#D95319', 'HandleVisibility', 'off');
    plot(time_vector, bs_s21, ':', 'LineWidth', 2.5, 'Color', '#EDB120', 'HandleVisibility', 'off');
    plot(time_vector, bs_s22, ':', 'LineWidth', 2.5, 'Color', '#7E2F8E', 'HandleVisibility', 'off');
    
    % Exact (Dashed Lines)
    plot(time_vector, ex_s11, '--', 'LineWidth', 1.5, 'Color', '#0072BD', 'HandleVisibility', 'off');
    plot(time_vector, ex_s12, '--', 'LineWidth', 1.5, 'Color', '#D95319', 'HandleVisibility', 'off');
    plot(time_vector, ex_s21, '--', 'LineWidth', 1.5, 'Color', '#EDB120', 'HandleVisibility', 'off');
    plot(time_vector, ex_s22, '--', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'HandleVisibility', 'off');
    
    plot(NaN, NaN, '-',  'LineWidth', 2,   'Color', [0.3 0.3 0.3],'DisplayName', 'GMRES');
    plot(NaN, NaN, ':',  'LineWidth', 2.5, 'Color', [0.3 0.3 0.3],'DisplayName', 'backslash');
    plot(NaN, NaN, '--', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3],'DisplayName', 'exact');

    grid on;
    xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Stress Value', 'FontSize', 12);
    title(sprintf('Stress Evolution at Point %d (F-AB-GMRES vs backslash vs Exact) (dt = %g)', target_pt, dt), ...
          'Interpreter', 'latex', 'FontSize', 15);
    legend('Location', 'best', 'NumColumns', 2); 
    
    % Subplot 2: Local Relative Error Evolution
    subplot(1, 2, 2);
    hold on;
    plot(time_vector, pt_rel_error_gmres, '-o', 'LineWidth', 2, 'Color', '#A2142F', 'MarkerFaceColor', '#A2142F', 'DisplayName', 'GMRES Error');
    plot(time_vector, pt_rel_error_bs, '-s', 'LineWidth', 2, 'Color', '#77AC30', 'MarkerFaceColor', '#77AC30', 'DisplayName', 'Backslash Error');
    grid on;
    xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', 12);
    
    % FIX 2: Use the exact matching matching LaTeX interpreter for the title to render norms nicely
    ylabel('Relative Error $\|\sigma_h - \sigma_{ex}\|_2 / \|\sigma_{ex}\|_2$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Local Relative Error at Point %d (dt = %g)', target_pt, dt), ...
          'Interpreter', 'latex', 'FontSize', 15);
    legend('Location', 'best');
end