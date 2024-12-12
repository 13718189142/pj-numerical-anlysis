function [eigenvalues, eigenvectors] = chebyshev_iteration(A, k, max_iter)
    n = size(A, 1);
    v = randn(n, 1);
    v = v / norm(v);
    for iter = 1:max_iter
        v = (A - 2 * mean(diag(A)) * eye(n)) * v;
        v = v / norm(v);
    end
    eigenvalues = v' * A * v;
    eigenvectors = v;
end
