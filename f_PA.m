% Parasitism probability
% PA is the parasitoid density
function f = f_PA(PA)

global k; % Interference parameter
global a; % Parasitoid attack rate

f = k * log(1.0 + a * PA / k);

end
