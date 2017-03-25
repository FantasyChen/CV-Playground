function F_matrices = calcFSevePoints(points1, points2)
% CALCULATE F MATRIX  USING SEVEN POINTS ALGORITHM
% RETURN MULTIPLE SOLUTIONS

A = [];
[T1, points1_norm] = normalization_2D(points1);
[T2, points2_norm] = normalization_2D(points2);
for i = 1:7
    point1 = points1_norm(:, i);
    point2 = points2_norm(:, i);
    A_i = kron(point2', point1');
    A = [A; A_i];
end

[~, ~, V] = svd(A);
F_1 = reshape(V(:,9), 3, 3)';
F_2 = reshape(V(:,8), 3, 3)';

syms x
result = double(vpa(solve(det(x*F_1 + F_2), 'Real', true)));

[m, ~] = size(result);
F_matrices = [];
count = 1;
for i = 1 : m
    if abs(imag(result(i))) < 0.0001
        F_matrices(:, :, count) = T2' * (real(result(i))*F_1 + F_2) * T1;
        count = count + 1;
    end
end

end