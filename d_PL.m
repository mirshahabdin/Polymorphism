% Parasitoid juvenile mortality
% Y parametrizes the degree of virulance of a parasitoid 
function d = d_PL(Y)

global dPL0; % Minimum parasitoid larval mortality
global cP; % cost of virulence

d = dPL0 * exp(cP * Y);

end