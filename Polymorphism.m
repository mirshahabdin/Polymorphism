globalConstants;

global N;
global M;

X = zeros(1, N);
Y = zeros(1, M);

% eta = eta_ij(X, Y);
eta = 0.0;
bHA = b_HA(X);
dPL = d_PL(Y);

global tauHA;
global tauHE;
global tauHL;
global tauHP;
global tauPA;
global tauPL;

delays = [tauHE, tauHL, tauHP, tauPL, tauHE+tauHL, tauHE+tauHL+tauHP, tauHE+tauHL+tauHP+tauHA, tauHA+tauHP, tauPL+tauPA];
hist = [];
poly = dde23('coevolution', delays, hist, [0, 5]);