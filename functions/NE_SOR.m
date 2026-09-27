function [z, res, p] = NE_SOR(A, v, omega, l)
%
%  Normal Equations Successive Over-Relaxation (NE-SOR) method
%  (Inner iteration routine)
%   this method applies a Kaczmarz-type relaxation sweep over the rows of A
%   to compute the preconditioned direction z approximating the solution to Az = v.
%   It is called internally during each outer loop of AB_GMRES.m.
%
%    INPUT:
%          A                : Matrix (Sparse format or not)
%          v                : right-hand side vector (current GMRES basis direction)
%          omega            : relaxation parameter
%          l                : number of allocation passes (inner iterations)
%
%    OUTPUT:
%          z                : preconditioned direction vector
%          res              : history of residual norms
%          p                : number of completed full sweeps
%
[m, n] = size(A);
z = zeros(n,1);
z_old = z;

res = norm(v - A*z);

if l == 0
    z = v; % return the starting direction
end

p = 0;

for p = 1:l
    for i = 1:m
        % Extract the i-th row of A and transpose it to a column vector
        alpha_i = A(i, :)';
        norm_alpha_sq = norm(alpha_i)^2;
        % Update the entire vector z
        z = z + omega * ((v(i) - alpha_i'*z) / norm_alpha_sq) * alpha_i;
    end
    if norm(z - z_old) < 1e-9
        % fprintf('Convergence with NE-SOR has been achieved!\n');
        return;
    end
    z_old = z;
    res = [res, norm(v - A*z)];
end