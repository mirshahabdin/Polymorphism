function v = test_func(t,y,Z)

ylag1 = Z(:,1);
ylag2 = Z(:,2)
v = zeros(3,1);
Z
v(1) = ylag1(1);
v(2) = ylag1(1) + ylag2(2);
v(3) = y(2);

end