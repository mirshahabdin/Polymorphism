function co = test_hist_coev(t, y, Z)

global N;

co = zeros(2*N, 1);

HL_D = zeros(N, 3);
HL = zeros(N, 1);
SHL_D = zeros(N, 3);
SHL = zeros(N, 1);

for i = 1:N
    HL_D(i, :) = Z(i, :);
    HL(i) = y(i);
    SHL_D(i, :) = Z(N+i, :);
    SHL(i) = y(N+i);
end

RHL = 0.0;
MHL = 0.0;
DHL = D_HL(HL, zeros(1, N));

co(1:N) = RHL - MHL - DHL;
co(N+1:2*N) = S_HL(zeros(N, 1), zeros(N, 1), HL, HL_D(:, 1), SHL);


end