% Matrix size
n = 500;

% Create a sparse, symmetric matrix
A = sprandn(n, n, density); % Sparse random matrix
A = (A + A') / 2; % Make it symmetric

% Number of iterations (can be adjusted for accuracy)
k = 50;
max_iter = 100;

% Measure time and compute eigenvalues for each method

% Lanczos
tic;
[eigenvalues_lanczos, eigenvectors_lanczos] = lanczos_algorithm(A, k);
time_lanczos = toc;

% Arnoldi
tic;
[eigenvalues_arnoldi, eigenvectors_arnoldi] = arnoldi_algorithm(A, k);
time_arnoldi = toc;

% Power Iteration
tic;
[eigenvalue_power, eigenvector_power] = power_iteration(A, max_iter);
time_power = toc;

% Inverse Iteration (using shift = 1)
tic;
[eigenvalue_inverse, eigenvector_inverse] = inverse_iteration(A, 1, max_iter);
time_inverse = toc;

% Subspace Iteration
tic;
[eigenvalues_subspace, eigenvectors_subspace] = subspace_iteration(A, k, 5);
time_subspace = toc;

% Chebyshev Iteration
tic;
[eigenvalues_chebyshev, eigenvectors_chebyshev] = chebyshev_iteration(A, k, max_iter);
time_chebyshev = toc;

% Display computation times and eigenvalues
fprintf('Lanczos method time: %.4f seconds\n', time_lanczos);
fprintf('Arnoldi method time: %.4f seconds\n', time_arnoldi);
fprintf('Power Iteration time: %.4f seconds\n', time_power);
fprintf('Inverse Iteration time: %.4f seconds\n', time_inverse);
fprintf('Subspace Iteration time: %.4f seconds\n', time_subspace);
fprintf('Chebyshev Iteration time: %.4f seconds\n', time_chebyshev);

% Display first few eigenvalues for comparison
fprintf('Lanczos eigenvalues: %f, %f, %f\n', eigenvalues_lanczos(1), eigenvalues_lanczos(2), eigenvalues_lanczos(3));
fprintf('Arnoldi eigenvalues: %f, %f, %f\n', eigenvalues_arnoldi(1), eigenvalues_arnoldi(2), eigenvalues_arnoldi(3));
fprintf('Power Iteration eigenvalue: %f\n', eigenvalue_power);
fprintf('Inverse Iteration eigenvalue: %f\n', eigenvalue_inverse);
 