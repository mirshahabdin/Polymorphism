% Larval host time-independent survival probabiliy
function sig = sigma_HL

global dHL; % Background host larval mortality
global tauHL; % Duration (in days) of host larval stage

sig = exp(-dHL * tauHL);

end