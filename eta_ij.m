% Probabiliy of encapsulation
% X parametrizes the degree of resistance in host
% Y parametrizes the level of virulence of the parasitoid
function eta = eta_ij(X, Y)

global A; % Slope of encapsulation trade-off curve
global M;
global N;

auxX = kron(X.', ones(1, M));
auxY = kron(ones(N, 1), Y);

eta = power(1 + exp(-A * (auxX - auxY)), -1);

end