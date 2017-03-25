function error = calcSampsonError(H, point1, point2)
% Calculation of the Sampson Error given two points correspondence and H
% ----------------------------------------------------------------------- %
H_vec = H';
H_vec = H_vec(:);
epsilon = zeros(2, 1);
epsilon = [zeros(1, 3), -point2(3) * point1', point2(2) * point1';
           point2(3)*point1', zeros(1, 3), -point2(1) * point1'] * H_vec;


J = zeros(2, 4);
J(1, 1) = -H_vec(4) + point2(2) * H_vec(7);
J(1, 2) = -H_vec(5) + point2(2) * H_vec(8);
J(1, 3) = 0;
J(1, 4) = point1(1) * H_vec(7) + point1(2) * H_vec(8) + H_vec(9);
J(2, 1) = H_vec(1) - point2(1) * H_vec(7);
J(2, 2) = H_vec(2) - point2(1) * H_vec(8);
J(2, 3) = -point1(1) * H_vec(7) - point1(2) * H_vec(8) - H_vec(9);

error = epsilon' * inv(J * (J')) * epsilon;

end