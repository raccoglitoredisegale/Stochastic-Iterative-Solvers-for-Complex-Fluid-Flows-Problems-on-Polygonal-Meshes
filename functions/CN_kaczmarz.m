function [results] = CN_kaczmarz(data_file, T, Nt, theta, choose_method, omega, tol, stagnation_tol, eta, outer_it_max, inner_it_max, Parameter_Tuning, adaptive)
%
%  Crank-Nicolson Time Integration of PDE system via F-AB-GMRES
%  (Time-stepping solver framework)
%
%   this method solves a time-dependent partial differential equation (PDE) system
%   using the Crank-Nicolson implicit scheme. At each temporal step, the resulting 
%   linear system is resolved via a Flexible AB-GMRES method preconditioned by 
%   Kaczmarz-type row projection inner solvers.
%
%    INPUT:
%          data_file        : string naming the .mat file containing system matrices
%          T                : total physical integration time
%          Nt               : total number of discrete time steps
%          theta            : parameter for time integration (0.5 for Crank-Nicolson)
%          choose_method    : - K for (relaxed) Kaczmarz
%                             - GK for Greedy Kaczmarz
%                             - RK for Randomized Kaczmarz
%                             - GRK for Greedy Randomized Kaczmarz
%          omega            : baseline relaxation parameter for inner solver
%          tol              : tolerance for the outer GMRES Krylov subspace method
%          stagnation_tol   : tolerance baseline to catch inner convergence stalls
%          eta              : tolerance passed down to inner Kaczmarz steps
%          outer_it_max     : maximum number of outer GMRES iterations per time step
%          inner_it_max     : maximum number of inner Kaczmarz row iterations
%          Parameter_Tuning : - 'Y' to re-tune parameters at every time step
%                             - 'F' to tune parameters only at the first time step
%                           
%          adaptive         : - 'Y' to utilize the adaptive restart variant
%                             - 'N' to run standard F-AB-GMRES configuration
%
%    OUTPUT:
%          sigma_sol        : solution vector at final time T across all nodes
%          sigma_history    : matrix tracking all solutions at every intermediate step
%          results          : structure collecting CPU/Wall time data, matrix conditions,
%                             residual profiles, and comparative quadrature errors
%

addpath('functions');
addpath('lymph_minimal')

fprintf('Loading data from %s...\n', data_file);
data = load(data_file);
results = struct();
results.adapt_struct = struct();

% extract variables from struct
M = data.OutputData.M;
A = data.OutputData.A;
Uold = data.OutputData.Uold;
F = data.OutputData.F;
femregion = data.OutputData.femregion;

% Time discretization parameters
dt  = T / Nt;
nts = floor(T / dt);
results.dt = dt;
results.T = T;

% Matrices for MX' + AX = F
ML_Fluid = M + A * theta * dt;
MR_Fluid = M - A * (1 - theta) * dt;
c = condition_number(ML_Fluid);
%results.c = condition_number((1/dt) * M + 0.5 * A);
results.c = c;
fprintf('Conditioning number with dt=%f is %e.      \n', dt, c);

% forcing vector at t=0
Fold = F.f_f * data.OutputData.source_sigma_t{1}(0) + F.g_f * data.OutputData.source_sigma_d_t{1}(0);

% this is where we will stor all the F_AB_GMRES solutions at every time step
sigma_history = zeros(length(Uold), nts + 1);
sigma_history(:, 1) = Uold;

% same but for backslash
Uold_bs = Uold;
backslash_history = zeros(length(Uold_bs), nts + 1);
backslash_history(:, 1) = Uold_bs;

% using LYMPH library structure
Solutions = cell(1, nts + 1);
Solutions{1} = SaveSolution(Uold, femregion);

results.time_wall_steps  = zeros(1, nts);
results.time_cpu_steps   = zeros(1, nts);
results.global_rel_error = zeros(1, nts + 1);

% initial error at t=0 using LYMPH structure
[Disp_init, ~] = GetSolutionQuadPointsPS(data.OutputData, femregion, Solutions{1}, 0);
num_quad_init = Disp_init.Sigma(:, 3:6); % column 3 to 6 are the F_AM_GMRES solutions
ex_quad_init  = Disp_init.Sigma(:, 7:10); % meanwhile 7 to 10 are the exact solutions

% tracking matrices
num_quad_points = size(Disp_init.Sigma, 1);
results.exact_quad_history = zeros(num_quad_points * 4, nts + 1);
results.num_quad_history   = zeros(num_quad_points * 4, nts + 1);
results.bs_quad_history    = zeros(num_quad_points * 4, nts + 1);
results.exact_quad_history(:, 1) = ex_quad_init(:);
results.num_quad_history(:, 1)   = num_quad_init(:);
results.bs_quad_history(:, 1)    = num_quad_init(:);
results.global_rel_error(1) = norm(num_quad_init(:) - ex_quad_init(:)) / norm(ex_quad_init(:));
results.global_rel_error_bs = zeros(1, nts + 1);

results.outer_iters = zeros(1, nts);
results.inner_iters = zeros(1, nts);
results.inner_avg = zeros(1,nts);

% timers
t = 0;
fprintf('\nStarting time: %g', t);
total_wall_time = tic;
total_cpu_time = cputime;
counter = 1;

% time loop
for t = dt : dt : nts*dt
    fprintf('\ntime: %g (Step %d/%d) using %s with: ', t, counter, nts, choose_method);
    Fnew = F.f_f * data.OutputData.source_sigma_t{1}(t) + F.g_f * data.OutputData.source_sigma_d_t{1}(t);

    % RHS for GMRES and Backslash
    F_Fluid    = MR_Fluid * Uold + dt * theta * Fnew + dt * (1 - theta) * Fold;
    F_Fluid_bs = MR_Fluid * Uold_bs + dt * theta * Fnew + dt * (1 - theta) * Fold;

    step_cpu_start = cputime;
    step_wall_start = tic;

    if Parameter_Tuning == 'Y'
        [inner_it_max, omega] = parameter_tuning(ML_Fluid, F_Fluid, eta, choose_method);
    elseif Parameter_Tuning == 'F'
        if t == dt
            [inner_it_max, omega] = parameter_tuning(ML_Fluid, F_Fluid, eta, choose_method);
        end
    end


    if adaptive == 'Y'
        [Uh, res_hist, adapt_struct] = adaptive_F_AB_GMRES(ML_Fluid, ...
            F_Fluid, Uold, choose_method, omega, tol, stagnation_tol, eta, outer_it_max, inner_it_max);
        results.adapt_struct = adapt_struct;
        results.inner_iters(counter) = adapt_struct.total_inner;
        results.outer_iters(counter) = adapt_struct.total_outer;
        inner_avg = adapt_struct.total_inner/adapt_struct.total_outer;
        fprintf('%d outer iterations, %d total inner iterations , %d inner_it/outer_it', adapt_struct.total_outer, adapt_struct.total_inner, inner_avg);

    else 
        [Uh, ~, res_hist, inner_iters, outer_iters] = F_AB_GMRES(ML_Fluid, ...
                F_Fluid, Uold, choose_method, omega, tol, eta, outer_it_max, inner_it_max);
        results.inner_iters(counter) = inner_iters;
        results.outer_iters(counter) = outer_iters;
        results.inner_avg(counter) = inner_iters/outer_iters;
        inner_avg = inner_iters/outer_iters;
        fprintf('%d outer iterations, %d total inner iterations , %d inner_it/outer_it', outer_iters, inner_iters, inner_avg);
    end

    results.res_history{counter} = res_hist;
    % BACKSLASH SOLVER
    Uh_bs = ML_Fluid \ F_Fluid_bs;
    results.time_wall_steps(counter) = toc(step_wall_start);
    results.time_cpu_steps(counter)  = cputime - step_cpu_start;

    % SaveSolution (from lymph) for GMRES
    SolStruct = SaveSolution(Uh, femregion);
    Solutions{counter + 1} = SolStruct;

    % SaveSolution (for Backslash)
    SolStruct_bs = SaveSolution(Uh_bs, femregion);

    % quadrature evaluation (also from lymph)
    [Disp, ~]    = GetSolutionQuadPointsPS(data.OutputData, femregion, {SolStruct}, t);
    [Disp_bs, ~] = GetSolutionQuadPointsPS(data.OutputData, femregion, {SolStruct_bs}, t);
    num_quad = Disp.Sigma(:, 3:6);
    bs_quad  = Disp_bs.Sigma(:, 3:6);
    ex_quad  = Disp.Sigma(:, 7:10);

    results.num_quad_history(:, counter + 1)   = num_quad(:);
    results.bs_quad_history(:, counter + 1)    = bs_quad(:);
    results.exact_quad_history(:, counter + 1) = ex_quad(:);

    % error Calculations
    norm_exact = norm(ex_quad(:), 2);
    norm_diff_gmres = norm(num_quad(:) - ex_quad(:), 2);
    norm_diff_bs    = norm(bs_quad(:) - ex_quad(:), 2);
    results.global_rel_error(counter + 1) = norm_diff_gmres / norm_exact;
    results.global_rel_error_bs(counter + 1) = norm_diff_bs / norm_exact;

    Uold    = Uh;
    Uold_bs = Uh_bs;
    Fold    = Fnew;

    sigma_history(:, counter + 1)     = Uh;
    backslash_history(:, counter + 1) = Uh_bs;

    counter = counter + 1;
end

results.time_wall_total = toc(total_wall_time);
results.time_cpu_total  = cputime - total_cpu_time;
sigma_sol = sigma_history(:, end);
results.sigma_sol = sigma_sol;
results.sigma_history = sigma_history;
results.Solutions = Solutions;
results.backslash_history = backslash_history; % Save raw backslash DoFs

fprintf(' Total Time : %.4f seconds (on average %.4f sec/step) \n', results.time_wall_total, results.time_wall_total/Nt);
end