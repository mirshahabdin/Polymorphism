% Adult host time-independent survival probabiliy
function sig = sigma_HA

global dHA; % Background host adualt mortality
global tauHA; % Duration (in days) of host adualt stage

sig = exp(-dHA * tauHA);

end