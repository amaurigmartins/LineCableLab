%% Quick example

A = [33 24; 48 57];

%% Correct way
[V,D] = eig(A);
A_sqrt_correct = V*sqrt(D)*inv(V)
E = V * diag(sqrt(diag(D))) / V
A_sqrt_correct*A_sqrt_correct

%% Wrong way
A_sqrt_wrong = sqrt(A)
A_sqrt_wrong*A_sqrt_wrong