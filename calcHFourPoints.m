function H = calcHFourPoints(points)
% Implementation of 4 points estimation of H Matrix. Here the input are
% homogeneous 2D points (3 x 1)
%-----------------------------------------------------------------------%
point1 = points(:, 1);
point2 = points(:, 2);
point3 = points(:, 3);
point4 = points(:, 4);

X = [point1, point2, point3];
Y = [point4];

lam = linsolve(X, Y);

H_inv = [lam(1) * point1, lam(2) * point2, lam(3) * point3];

H = inv(H_inv);


end