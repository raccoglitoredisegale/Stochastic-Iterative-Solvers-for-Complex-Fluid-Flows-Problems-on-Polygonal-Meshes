function cond_number = condition_number(A)
%
%   this method calculates the 2-norm condition number of a matrix A 
%    INPUT:
%          A                : Matrix (Sparse format or not)
%
%    OUTPUT:
%          cond_number      : 2-norm condition number of the matrix
%
S = svd(full(A), 'econ');

sigma_max = max(S);
sigma_min = min(S);

cond_number = sigma_max / sigma_min;
end
