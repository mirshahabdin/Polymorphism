% Adult parasitiod death rate
% HL density of host larvae
% PA is the parasitoid density
function R = R_PA(HL, PA)

global M;
global eta;
global eps;

R = zeros(1, M);
aux = (HL * (1.0 - eta)) .* f_PA(PA) * sigma_PL;

if M > 1
    for i = 1:M
        R(i) = (1.0 - eps) * aux(i) + eps / (M - 1) * sum(aux(setdif(1:end, i)));
    end
else
    R = aux;
end

end