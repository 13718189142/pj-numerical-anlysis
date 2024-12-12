function [eigenvalue, eigenvector] = inverse_iteration(A, shift, max_iter)
    n = size(A, 1);
    I = speye(n);
    b = randn(n, 1);
    b = b / norm(b);
    for i = 1:max_iter
        b = (A - shift * I) \ b;
        b = b / norm(b);
    end
    eigenvalue = b' * A * b;
    eigenvector = b;
end
