function h = history(t, alpha, gamma, dA, dP)

h = zeros(4, 1);

global P0;
global A0;
global tA;

h(1) = gamma*A0/(alpha*P0)*(1-exp(-alpha*P0*tA));
h(2) = A0;
h(3) = P0;
h(4) = tA*P0;