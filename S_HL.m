% Larval host time-dependent survival probability
function S = S_HL(PAT, PAt, HLT, HLt, SHL)

global eta;
global c;

S = -((1.0 - eta) * (f_PA(PAT) - f_PA(PAt)).' + c * (HLT - HLt)) .* SHL;

end