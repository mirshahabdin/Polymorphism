% Adult parasitoid death rate
% PA is the parasitoid density
function D = D_PA(PA)

global dPA; % Background parasitoid adult mortality

D = dPA * PA;

end