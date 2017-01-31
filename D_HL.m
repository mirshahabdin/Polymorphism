% Larval host death rate
% HL density of host larvae
% PA is the parasitoid density
function D = D_HL(HL, PA)

%global c; % Competition mortality coefficient\
c = 0.0;
global dHL; % Background host larval mortality
global eta;

D = ((1.0 - eta) * (f_PA(PA).') + c * HL + dHL) .* HL;

end