function [R, t] = umeyama(points_3D, points_3D_cam)
% UMEYAMA_METHOD
% Implementation of Umeyama's method

X = points_3D;
Y = points_3D_cam;

mu_X = mean(X, 2);
mu_Y = mean(Y, 2);
   

var_X = 0;
var_Y = 0;
for j = 1:size(X, 2)
    var_X = var_X + (norm(X(:,j)-mu_X))^2;
    var_Y = var_Y + (norm(Y(:,j)-mu_Y))^2;
end

var_X = var_X /size(X, 2);
var_Y = var_Y /size(X, 2);

covXY = zeros(3, 3);
for j = 1:size(X, 2)
    covXY = covXY + (Y(:, j) - mu_Y) * (X(:, j) - mu_X)';
end
covXY = covXY/size(X, 2);
[U, D, V] = svd(covXY);
% need to check rank of covXY
S = eye(3);
if(det(U) * det(V) == -1)
    S(3, 3) = -1;
end
R = U * S * V';
%c = 1/var_X * trace(D * S);
t = mu_Y - R * mu_X;


end