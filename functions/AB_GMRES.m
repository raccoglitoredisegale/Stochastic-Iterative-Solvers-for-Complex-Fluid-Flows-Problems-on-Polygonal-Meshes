function [u,res,total_inner,avg_inner]=AB_GMRES(A,b,u_0,omega,tol,m_max,l_max)
%
%  AB-GMRES method preconditioned by inner iterations
%  (algorithm 2.1)
%   this method solves a linear system Au=b using GMRES
%   method starting at the initial guess u0 up to a tolerance tol using
%   at most m_max iterations.
%
%    INPUT:
%          A                : Matrix (Sparse format or not)
%          b                : vector
%          u_0              : initial solution
%          omega            : parameter for NE-SOR (if omitted, default=1)
%          tol              : tolerance (if omitted, default is 1e-6)
%          m_max            : maximum number of outer iterations
%                             (if omitted, m_max = 100)
%          l_max            : maximum number of inner NE-SOR iterations
%                             (if omitted, l_max = 15)
%
%    OUTPUT:
%          u                : solution vector
%          res              : history of norms of the residuals
%          total_inner      : total number of inner row updates performed
%          avg_inner        : average number of inner iterations per outer loop
%

if nargin < 7, l_max = 15; end            % default values
if nargin < 6, m_max = 100; end
if nargin < 5, tol = 1e-6; end
if nargin < 4, omega = 1; end

r_0 = b - A*u_0;
beta = norm(r_0);
res = beta;                % Array to store residual history
V(:,1) = r_0 / beta;        % Matrix to store the searching directions (v_1)
H = [];                    % Storing the Hessenberg matrix

k = 0;                     % number of (outer) iterations
norm_b = norm(b);

m = size(A, 1);
total_inner = 0; % Track total row updates

% outer iterations loop
while (res(end) / norm_b > tol) && (k < m_max)
    k = k + 1;

    % Aplly NE-SOR algorithm
    [z_k, ~, p_k] = NE_SOR(A, V(:,k), omega, l_max);
    %z_k = V(:,k);
    total_inner = total_inner + (p_k * m);
    w_k = A*z_k;



    % Orthogonalize the new direction w_k against all previous directions
    for i = 1:k
        H(i,k) = w_k' * V(:,i);
        w_k = w_k - H(i,k) * V(:,i);
    end
    H(k+1,k) = norm(w_k);
    V(:,k+1) = w_k / H(k+1,k); % am i sure it is not zero?

    % Check if the matrix is singular. if it is, the algorithm has stalled
    % the new directions are not linearly independent anymore. should we
    % break or restart? or what?
    if rank(H(1:k+1, 1:k)) < k
        warning('NE-SOR: H matrix has become singular.');
        %then we save the last solution at k-1.
        e_1 = zeros(k,1);
        e_1(1) = beta;
        y_k = H(1:k, 1:(k-1)) \ e_1;
        u_k = V(:, 1:(k-1)) * y_k;
        y = NE_SOR(A, u_k, omega, l_max);
        z = A' * y;
        u = u_0 + z;
        return;
    end
    % Minimize ||beta*e1 - H(1:k+1, 1:k)*y||
    % (metodo minimi quadrati)
    e_1 = zeros(k+1,1);
    e_1(1) = beta;
    y_k = H(1:k+1, 1:k) \ e_1;
    %[Q, R] = qr(H(1:k+1, 1:k), 0);
    %y_k = R \ (Q' * e_1);

    % Update and store residual
    res(k + 1) = norm(e_1 - H(1:k+1, 1:k) * y_k);

end

% u_k = [v_1, ..., v_k] * y_k
u_k = V(:, 1:k) * y_k;

[z, ~, p_final] = NE_SOR(A, u_k, omega, l_max);
total_inner = total_inner + (p_final * m);
u = u_0 + z;
avg_inner = round(total_inner / k);

% fprintf(' NE-SOR:  %d iterations completed \n',k );


end
