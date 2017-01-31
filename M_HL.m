% Larval host recruitment rate
% HA is density of host larvae
function M = M_HL(RHL, SHL)

M = RHL .* SHL * sigma_HL;

end