% Host fecundity
% X parametrizes the degree of resistance in host
function b = b_HA(X)

global b0; % Maximum host fecunidty
global cH; % Cost of resistance

b = b0 * exp(-cH * X);

end