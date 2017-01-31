% Adult parasitoid time-independent survival probabiliy
function sig = sigma_PA

global dPA; % Background parasitoid adult mortality
global tauPA; % Duration (in days) of parasitoid adult stage

sig = exp(-dPA * tauPA);

end