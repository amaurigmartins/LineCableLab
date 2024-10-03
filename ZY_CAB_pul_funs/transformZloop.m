% transformZloop.m
function Z_prime = transformZloop(Z)
% Transforms loop-based impedances into conductor quantities
% (core/sheath/armor)

% Check the size of the Z matrix
N = size(Z, 1); % N is the number of layers, assuming Z is NxN

% Build the voltage transformation matrix T_V
T_V = eye(N, N+1); % Start with an identity matrix
for i = 1:N
    T_V(i, i+1) = -1; % Set the -1 in the next column
end
T_V(:,N+1)=[];

% Build the current transformation matrix T_I
T_I = tril(ones(N)); % Lower triangular matrix of ones

% Compute the new impedance matrix Z_prime
Z_prime = T_V \ Z * T_I;
end