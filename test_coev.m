function co = test_coev(t, y, Z)

global N;

co = zeros(3*N, 1);

HL_D = zeros(N, 3);
HL = zeros(N, 1);
HA_D = zeros(N, 3);
HA = zeros(N, 1);
SHL_D = zeros(N, 3);
SHL = zeros(N, 1);

for i = 1:N
    HL_D(i, :) = Z(i, :);
    HL(i) = y(i);
    HA_D(i, :) = Z(N+i, :);
    HA(i) = y(N+i);
    SHL_D(i, :) = Z(2*N+i, :);
    SHL(i) = y(2*N+i);
end

RHL = R_HL(HA);
MHL = M_HL(R_HL(HA_D(:, 1)), SHL);
DHL = D_HL(HL, zeros(1, N));

RHA = R_HA(MHL); 
MHA = zeros(1, N);
DHA = D_HA(HA_D(:, 3), SHL_D(:, 2));

co(1:N) = RHL - MHL - DHL;
co(N+1:2*N) = RHA - MHA - DHA;
co(2*N+1:3*N) = S_HL(zeros(N, 1), zeros(N, 1), HL, HL_D(:, 1), SHL);

% if (t == 28)
%     co(N+1:2*N) = 15.0 * SHL;
% end

end