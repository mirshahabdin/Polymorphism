function co = test(t, y, Z)

global bHA;
global dHL;
global dHA;
global c;

co = zeros(3, 1);


HL_D = Z(1, :);
HL = y(1);
HA_D = Z(2, :);
HA = y(2);
SHL_D = Z(3, :);
SHL = y(3);

co(1) = bHA * HA_D(1) * sigma_HE - bHA * HA_D(2) * sigma_HE * SHL * sigma_HL - dHL * HL;
co(2) = bHA * HA_D(4) * sigma_HE * SHL_D(3) * sigma_HL * sigma_HP - ...
    bHA * HA_D(6) * sigma_HE * SHL_D(5) * sigma_HL * sigma_HP * sigma_HA - dHA * HA;
co(3) = - c * (HL - HL_D(7)) * SHL;


end