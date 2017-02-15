globalConstants;

global eta;
global tauHL;
global tauHE;
global tauHA;
global tauHP;

eta = 1;

delays = [tauHE, tauHE+tauHL, tauHP, tauHE+tauHL+tauHP, tauHA+tauHP, tauHE+tauHL+tauHP+tauHA, tauHL];

poly = dde23('test', delays, 'test_hist', [0, 5000]);