function points_3D_inliers_proj_dehomo = ePNP(points_2D_inliers, points_3D_inliers, K)

num_inliers = size(points_3D_inliers, 2);
points_2D_inliers_norm = inv(K) * [points_2D_inliers; ones(1, num_inliers)];

mean_3D = mean(points_3D_inliers, 2);
cov_3D = (points_3D_inliers - mean_3D)*(points_3D_inliers - mean_3D)'/(num_inliers - 1);

[U, D, V] = svd(cov_3D);
d1 = D(1, 1);
d2 = D(2, 2);
d3 = D(3, 3);
s = sqrt((d1+d2+d3)/3);
v1 = V(:, 1);
v2 = V(:, 2);
v3 = V(:, 3);

% Define control points
c1 = mean_3D;
c2 = s * v1 + mean_3D;
c3 = s * v2 + mean_3D;
c4 = s * v3 + mean_3D;
A = [c2 - c1, c3 - c1, c4 - c1];
points_3D_inliers_param = zeros(size(points_3D_inliers));
alpha_matrix = zeros(4, num_inliers); 
% Parameterization of 3D points
for i = 1:num_inliers
    cur_point = points_3D_inliers(:, i);
    alpha = inv(A) * (cur_point - c1);
    alpha_matrix(2, i) = alpha(1);
    alpha_matrix(3, i) = alpha(2);
    alpha_matrix(4, i) = alpha(3);
    alpha_matrix(1, i) = 1 - alpha(1) - alpha(2) - alpha(3);
    points_3D_inliers_param(:, i) = alpha_matrix(1, i) * c1 + alpha_matrix(2, i) * c2 + alpha_matrix(3, i) * c3 + alpha_matrix(4, i) * c4;
    points_3D_inliers_param(:, i) = points_3D_inliers_param(:, i) / points_3D_inliers_param(3, i);
end

% Define matrix M
M = zeros(2 * num_inliers, 12);
for i = 1:num_inliers
    M_cur = [];
    for j = 1:4
        temp = [alpha_matrix(j, i), 0, -alpha_matrix(j, i) * points_2D_inliers_norm(1, i);
                0, alpha_matrix(j, i), -alpha_matrix(j, i) * points_2D_inliers_norm(2, i)];
        M_cur = [M_cur, temp];
    end
    M = [M; M_cur];
end

[~, ~, V] = svd(M);
control_global_vec = V(:, 12);
control1 = control_global_vec(1:3);
control2 = control_global_vec(4:6);
control3 = control_global_vec(7:9);
control4 = control_global_vec(10:12);

points_3D_inliers_cam = zeros(3, num_inliers);
for i = 1: num_inliers
    points_3D_inliers_cam(:, i) = alpha_matrix(1, i) * control1 + alpha_matrix(2, i) * control2 + ...
                        alpha_matrix(3, i) * control3 + alpha_matrix(4, i) * control4;
end
mean_3D_cam = mean(points_3D_inliers_cam, 2);
sigma_3D_cam = var(points_3D_inliers_cam,0,2);
sigma_3D_cam = sum(sigma_3D_cam);
sigma_3D = var(points_3D_inliers, 0, 2);
sigma_3D = sum(sigma_3D);


for i = 1: num_inliers
    if points_3D_inliers_cam(3, i) < 0
        points_3D_inliers_cam(:, i) = - points_3D_inliers_cam(:, i) * sqrt(sigma_3D/sigma_3D_cam);
    else
        points_3D_inliers_cam(:, i) = points_3D_inliers_cam(:, i) * sqrt(sigma_3D/sigma_3D_cam);
    end
end

[R, t] = umeyama(points_3D_inliers, points_3D_inliers_cam);

R = min_model(1:3, 1:3);
t =  min_model(1:3, 4);
points_3D_inliers_proj = K * [R, t] * [points_3D_inliers; ones(1, num_inliers)];
for i =1:num_inliers
    points_3D_inliers_proj_dehomo(:, i) = points_3D_inliers_proj(:, i)/points_3D_inliers_proj(3, i);
end
