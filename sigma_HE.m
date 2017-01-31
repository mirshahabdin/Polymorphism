% Egg host time-independent survival probabiliy
function sig = sigma_HE

global dHE; % Background host egg mortality
global tauHE; % Duration (in days) of host egg stage

sig = exp(-dHE * tauHE);

end