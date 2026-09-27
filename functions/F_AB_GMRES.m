function [u, u_history, res, total_inner, total_outer, avg_inner] = F_AB_GMRES(A,b,u_0,choose_method,omega,tol,eta,m_max,l_max)
%
%  Flexible AB-GMRES preconditioned by Kaczmarz-type methods
%  (algorithm 3.1)
%
%   [u,u_history,res,total_inner,total_outer,avg_inner]=F_AB_GMRES(A,b,u_0,choose_method,omega,tol,eta,m_max,l_max)
%
%   this method solves a linear system Au=b using a flexible AB-preconditioned GMRES
%   method starting at the initial guess u0 up to a tolerance tol using
%   at most m_max iterations. the function will call kaczmark_solve.m to
%   apply the desired algorithm for resolving every inner iteration.
%
%    INPUT:
%          A                : Matrix (Sparse format or not)
%          b                : vector
%          u_0              : initial solution
%          choose_method    : - K for (relaxed) Kaczmarz
%                             - GK for Greedy Kaczmarz
%                             - RK for Randomized Kaczmarz
%                             - GRK for Greedy Randomized Kaczmarz
%          omega            : parameter for Kaczmarz function
%                             (if omitted, default is 1)
%          tol              : tolerance for outer loop stagnation
%                             (if omitted, default is 1e-6)
%          eta              : tolerance for inner loop solver 
%                             (if omitted, default is 1e-6)
%          m_max            : maximum number of outer iterations
%                             (if omitted, m_max = 100)
%          l_max            : maximum number of inner iterations of the
%                             selected K. method (if omitted, l_max = 100)
%
%    OUTPUT:
%          u                : solution vector
%          u_history        : matrix of iterates (solutions at each k)
%          res              : history of norms of the residuals
%          total_inner      : total number of inner iterations
%          total_outer      : total number of outer iterations
%          avg_inner        : average number of inner iterations per outer iteration
%

if nargin < 9, l_max = 100; end            % default values
if nargin < 8, m_max = 100; end
if nargin < 7, eta = 1e-6; end
if nargin < 6, tol = 1e-6; end
if nargin < 5, omega = 1; end

r_0 = b - A*u_0;
beta = norm(r_0);
res = beta;                % Array to store residual history
V(:,1) = r_0 / beta;        % Matrix to store the searching directions
H = [];                    % Storing the Hessenberg matrix
Z = [];                   % Where we store the preconditioned directions
u_history(:,1) = u_0;           % Where we store solution vector at each iter
k = 0;                     % number of (outer) iterations
norm_b = norm(b);
total_inner = 0;

% outer iterations loop
while (res(end) / norm_b > tol) && (k < m_max)
    k = k + 1;
    % fprintf('outer iteration: %d \n', k);

    % inner iterations loop: l_k iterations of a Kaczmark-type method
    % to solve A*z_k = V(:,k) to obtain z_k.
    [z_k, ~, l_k, ~] = kaczmarz_solve (A, V(:,k), l_max, eta, choose_method, omega);
    total_inner = total_inner + l_k;

    % z_k = A' * V(:,k);
    Z(:,k) = z_k;
    w_k = A*z_k;

    % Orthogonalize the new direction w_k against all previous directions
    for i = 1:k
        H(i,k) = w_k' * V(:,i);
        w_k = w_k - H(i,k) * V(:,i);
    end
    H(k+1,k) = norm(w_k);

    if rcond(H(1:k, 1:k)) < 1e-12
        fprintf('Fatal breakdown: H became singular at iteration %d using %s. Stagnation reached.\n', k, choose_method);
        % Return the LAST valid solution before the math blew up
        u = u_history(:, end);
        total_outer = k;
        avg_inner = round(total_inner/k);
        return;
    end

    if H(k+1,k) < 1e-9
        fprintf('Found exact sol.  in %d outer it \n', k);
        % Solve the least squares problem one last time
        e_1 = zeros(k, 1); % Notice we only need size k now
        e_1(1) = beta;
        % We solve using the square k x k matrix (H_k)
        y_k = H(1:k, 1:k) \ e_1;
        u_k = u_0 + Z(:, 1:k) * y_k;
        u_history = [u_history, u_k];
        res(k + 1) = 0; % Residual is effectively zero
        u = u_history(:, end);
        total_outer = k;
        avg_inner = round(total_inner/k);
        return;
    end

    % Standard Least Squares Solve (for when we haven't broken down)
    e_1 = zeros(k+1, 1);
    e_1(1) = beta;
    y_k = H(1:k+1, 1:k) \ e_1;

    % Update and append current correction
    u_k = u_0 + Z(:, 1:k) * y_k;
    u_history = [u_history, u_k];

    % Update and store residual
    res(k + 1) = norm(e_1 - H(1:k+1, 1:k) * y_k);

    % Prepare next basis vector
    V(:,k+1) = w_k / H(k+1,k);
end

total_outer = k;
% fprintf("Loop of %s finished in %d outer iter\n", choose_method, k);
u = u_history(:, end);
avg_inner = round(total_inner/k);
end