%delays = [tauHE, tauHL, tauHP, tauPL, tauHE+tauHL, tauHE+tauHL+tauHP, tauHE+tauHL+tauHP+tauHA, tauHA+tauHP, tauPL+tauPA];
function co = coevolution(t, y, Z)

global N;
global M;

co = zeros(3*N+M, 1);

HL_D = zeros(9, N);
HL = zeros(1, N);
HA_D = zeros(9, N);
HA = zeros(1, N);
SHL_D = zeros(9, N);
SHL = zeros(1, N);
PA_D = zeros(9, M);
PA = zeros(1, M);

for i = 1:N
    HL_D(:, i) = Z(:, i);
    HL(i) = y(i);
    HA_D(:, i) = Z(:, N+i);
    HA(i) = y(N+i);
    SHL_D(:, i) = Z(:, 2*N+i);
    SHL(i) = y(2*N+i);
end

for i = 1:M
    PA_D(:, i) = Z(:, 3*N+i);
    PA(i) = y(3*N+i);
end


RHL = R_HL(HA_D(1, :));
MHL = M_HL(R_HL(HA_D(5, :)), SHL);
DHL = D_HL(HL, PA);

RHA = R_HA(M_HL(R_HL(HA_D(6, :)))); 
MHA = M_HA(R_HA(M_HL(R_HL(HA(7, :)),SHL_D(8, :))));
DHA = D_HA(HA);

RPA = R_PA(HL_D(4, :), PA_D(4, :));
MPA = M_PA(R_PA(HL_D(9, :), PA_D(9, :)));
DPA = D_PA(PA);


co(1:N) = HL(RHL, MHL, DHL);
co(N+1:2*N) = HA(RHA, MHA, DHA);
co(2*N+1:3*N) = S_HL(PA, PA_D(2, :), HL, HL_D(2, :));
co(3*N:end) = PA(RPA, MPA, DPA);


end