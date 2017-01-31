% Adult host death rate
% HA is density of host larvae
function D = D_HA(HA)
%function D = D_HA(HA, SHL)

global dHA; % Background host adualt mortality
%global bHA;

D = dHA * HA;
%D = bHA * HA .* SHL;

end