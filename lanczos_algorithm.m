function [eigenvalues, eigenvectors] = lanczos_algorithm(A, k)
    % Input:
    % A - The large sparse symmetric matrix (size n x n)
    % k - The number of Lanczos iterations (usually much smaller than n)
    %
    % Output:
    % eigenvalues - The computed eigenvalues of the tridiagonal matrix
    % eigenvectors - The computed eigenvectors

    % Get matrix size
    n = size(A, 1);

    % Step 1: Initialization
    % Random vector v_1
    v = randn(n, 1);
    v = v / norm(v);  % Normalize the vector

    % Initialize alpha and beta
    alpha = zeros(k, 1);
    beta = zeros(k, 1);
    V = zeros(n, k);  % Store the orthonormal vectors
    T = zeros(k, k);  % The tridiagonal matrix

    V(:,1) = v;

    % Step 2: Lanczos Iterations
    for j = 1:k
        % Step 2.1: Compute w_j = A * v_j
        w = A * V(:,j);

        % Step 2.2: Orthogonalize w_j with respect to v_(j-1)
        if j > 1
            w = w - beta(j-1) * V(:,j-1);
        end

        % Step 2.3: Compute alpha_j = v_j^T * w_j
        alpha(j) = V(:,j)' * w;

        % Step 2.4: Compute the new vector w_j
        w = w - alpha(j) * V(:,j);

        % Step 2.5: Compute beta_j = norm(w_j)
        beta(j) = norm(w);

        % Step 2.6: If beta_j is small, stop early (break)
        if beta(j) < 1e-10
            break;
        end

        % Step 2.7: Normalize the vector to get the next v_j
        V(:,j+1) = w / beta(j);
    end

    % Step 3: Construct the tridiagonal matrix T
    T(1:k, 1:k) = diag(alpha) + diag(beta(2:k), 1) + diag(beta(2:k), -1);

    % Step 4: Compute the eigenvalues of the tridiagonal matrix T
    [Q, D] = eig(T);
    eigenvalues = diag(D);  % Extract eigenvalues
    eigenvectors = V(:,1:k) * Q;  % Eigenvectors of A

end

