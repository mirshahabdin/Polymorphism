% Larval host recruitment rate
% HA is density of host larvae
function R = R_HL(HA)

global bHA; % Daily host adult fecundity (number of eggs)
global eps;
global N;

R = zeros(1, N);
aux = bHA .* HA * sigma_HE;

if N > 1
    for i = 1:N
        R(i) = (1 - eps) * aux(i) + eps / (N - 1) * sum(aux(setdiff(1:end, i)));
    end
else
    R = aux;
end

end