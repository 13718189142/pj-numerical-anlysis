function [eigenvalue, eigenvector] = power_iteration(A, max_iter)
    n = size(A, 1);
    b = randn(n, 1);
    b = b / norm(b);
    for i = 1:max_iter
        b = A * b;
        b = b / norm(b);
    end
    eigenvalue = b' * A * b;
    eigenvector = b;
end
