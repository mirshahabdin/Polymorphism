global P0;
global A0;
global tA;

P0 = 0.3;
A0 = 1.45;
tA = 1;
tP = 1;
alpha = 1;
gamma = 1.39;
dA = 0.4;
dP = 0.55;

sol = dde23('mathematica', [tA, tP], 'history', [0, 50], [], alpha, gamma, dA, dP);

figure;
plot(sol.x, sol.y(1:3, :));