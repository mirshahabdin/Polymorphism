function co = mathematica(t, y, Z, alpha, gamma, dA, dP)

co = zeros(4, 1);

colag1 = Z(:, 1);
colag2 = Z(:, 2);

co(1) = gamma*y(2) - gamma*colag1(2)*exp(-alpha*y(4)) - alpha*y(3)*y(1);
co(2) = gamma*colag1(2)*exp(-alpha*y(4)) - dA*y(2);
co(3) = alpha*colag2(3)*colag2(1) - dP*y(3);
co(4) = y(3)-colag1(3);

end
