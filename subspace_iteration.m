function [eigenvalues, eigenvectors] = subspace_iteration(A, k, num_eigenvalues)
    n = size(A, 1);
    V = randn(n, k);
    for iter = 1:100
        B = A * V;
        [Q, R] = qr(B, 0);
        V = Q;
        eigenvalues = eig(R);
        if length(unique(eigenvalues)) == num_eigenvalues
            break;
        end
    end
    eigenvectors = V;
end
