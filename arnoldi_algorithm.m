function [eigenvalues, eigenvectors] = arnoldi_algorithm(A, k)
    n = size(A, 1);
    v = randn(n, 1);
    v = v / norm(v);
    H = zeros(k, k);
    V = zeros(n, k);
    V(:, 1) = v;
    
    for j = 1:k
        w = A * V(:, j);
        for i = 1:j
            H(i, j) = V(:, i)' * w;
            w = w - H(i, j) * V(:, i);
        end
        H(j + 1, j) = norm(w);
        if H(j + 1, j) > 1e-10
            V(:, j + 1) = w / H(j + 1, j);
        else
            break;
        end
    end
    [Q, D] = eig(H(1:k, 1:k));
    eigenvalues = diag(D);
    eigenvectors = V(:, 1:k) * Q;
end
