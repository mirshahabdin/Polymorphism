% Pupal host time-independent survival probabiliy
function sig = sigma_HP

global dHP; % Background host pupal mortality
global tauHP; % Duration (in days) of host pupal stage

sig = exp(-dHP * tauHP);

end