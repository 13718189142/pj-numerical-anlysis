% Parameters
n = 1000; % Size of the matrix
density = 0.01; % Fraction of non-zero entries

% Generate a random sparse symmetric matrix
A = sprandn(n, n, density); % Sparse random matrix
A = (A + A') / 2; % Make it symmetric

% Compute the eigenvalues
eigenvalues = eigs(A, 20); % Compute 20 largest magnitude eigenvalues

% Plot eigenvalues
figure;
histogram(eigenvalues, 30);
title('Eigenvalue Distribution');
xlabel('Eigenvalue');
ylabel('Frequency');
