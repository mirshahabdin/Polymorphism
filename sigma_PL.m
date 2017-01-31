% Juvenile parasitoid time-independent survival probabiliy
function sig = sigma_PL

global dPL; % Background parasitoid juvenile mortality
global tauPL; % Duration (in days) of parasitoid egg and larval stage

sig = exp(-dPL * tauPL);

end