function [points_3D_homo, P_prime] = my_triangulate(points1, points2, F)
num_points = size(points1, 2);
points_3D_homo = zeros(4, num_points);


[U,D,V] = svd(F);
D_prime = D;
D_prime(3,3) = (D(1,1)+D(2,2))/2;

Z = zeros(3,3);
Z(1,2) = -1;
Z(2,1) = 1;
Z(3,3) = 1;

W = zeros(3,3);
W(1,2) = 1;
W(2,1) = -1;

M = U*Z*D_prime*(V');
S = U*W*U';
e(1) = S(3,2);
e(2) = S(1,3);
e(3) = S(2,1);


P_prime = zeros(3,4);
P_prime(:,1:3) = M;
P_prime(:,4) = e;


for i = 1 : num_points
    cur_point1 = points1(:, i);
    cur_point2 = points2(:, i);
    
    T = zeros(3,3);
    T(1,1) = cur_point1(3);
    T(2,2) = cur_point1(3);
    T(3,3) = cur_point1(3);
    T(1,3) = (-1)*cur_point1(1);
    T(2,3) = (-1)*cur_point1(2);
    
    T_prime = zeros(3,3);
    T_prime(1,1) = cur_point2(3);
    T_prime(2,2) = cur_point2(3);
    T_prime(3,3) = cur_point2(3);
    T_prime(1,3) = (-1)*cur_point2(1);
    T_prime(2,3) = (-1)*cur_point2(2);
    
    F_s = (inv(T_prime'))*F*inv(T);
    
    [~,~,V] = svd(F_s);
    e_s = V(:,3);
    [~,~,V] = svd(F_s');
    e_s_prime = V(:,3);
    
    e_s = e_s*sqrt(1/(e_s(1)^2 + e_s(2)^2));
    e_s_prime = e_s_prime*sqrt(1/(e_s_prime(1)^2 + e_s_prime(2)^2));
    
    R = zeros(3,3);
    R(1,1) = e_s(1);
    R(1,2) = e_s(2);
    R(2,1) = (-1)*e_s(2);
    R(2,2) = e_s(1);
    R(3,3) = 1;
    
    R_prime = zeros(3,3);
    R_prime(1,1) = e_s_prime(1);
    R_prime(1,2) = e_s_prime(2);
    R_prime(2,1) = (-1)*e_s_prime(2);
    R_prime(2,2) = e_s_prime(1);
    R_prime(3,3) = 1;
    
    F_s = R_prime*F_s*R';
    
    a = F_s(2,2);
    b = F_s(2,3);
    c = F_s(3,2);
    d = F_s(3,3);
    f = e_s(3);
    f_prime = e_s_prime(3);
    
    syms t
    s = double(vpa(solve(t*((a*t+b)^2 + f_prime*((c*t+d)^2))^2 - (a*d - b*c)*((1 + (f*t)^2)^2)*(a*t+b)*(c*t+d), 'Real', true)));
    
    best_value = Inf;
    best_t = Inf;
    for j = 1:size(s)
        t = real(s(j));
        cur_value = (t^2)/(1 + (f*t)^2) + ((c*t+d)^2)/((a*t+b)^2 + (f_prime^2)*((c*t+d)^2));
        if(cur_value < best_value)
            best_t = t;
            best_value = cur_value;
        end
    end
    
    t = best_t;
    
    l_vec = zeros(3,1);
    l_vec(1) = t*f;
    l_vec(2) = 1;
    l_vec(3) = (-1)*t;
    
    l_vec_prime = zeros(3,1);
    l_vec_prime(1) = (-1)*f_prime*(c*t+d);
    l_vec_prime(2) = a*t+b;
    l_vec_prime(3) = c*t+d;
    
    cur_point1_hat = zeros(3,1);
    cur_point1_hat(1) = (-1)*l_vec(1)*l_vec(3);
    cur_point1_hat(2) = (-1)*l_vec(2)*l_vec(3);
    cur_point1_hat(3) = l_vec(1)^2 + l_vec(2)^2;
    
    
    cur_point2_hat = zeros(3,1);
    cur_point2_hat(1) = (-1)*l_vec_prime(1)*l_vec_prime(3);
    cur_point2_hat(2) = (-1)*l_vec_prime(2)*l_vec_prime(3);
    cur_point2_hat(3) = l_vec_prime(1)^2 + l_vec_prime(2)^2;
    
    cur_point1_hat = inv(T)*(R')*cur_point1_hat;
    cur_point1_hat = cur_point1_hat/cur_point1_hat(3);
    
    cur_point2_hat = inv(T_prime)*(R_prime')*cur_point2_hat;
    cur_point2_hat = cur_point2_hat/cur_point2_hat(3);

    l_prime_vec = F * cur_point1_hat;
    l_prime_orth = zeros(3,1);
    l_prime_orth(1) = (-1)*l_prime_vec(2)*cur_point2_hat(3);
    l_prime_orth(2) = l_prime_vec(1)*cur_point2_hat(3);
    l_prime_orth(3) = l_prime_vec(2)*cur_point2_hat(1) - l_prime_vec(1)*cur_point2_hat(2);
    
    P_i_vec = (P_prime')*l_prime_orth;
    
    P_eye = zeros(3,4);
    P_eye(:,1:3) = eye(3);
    
    [~, ~, V] = svd(P_eye);
    
    point_3D_1 = V(:,4);

    point_3D_2 = (P_eye')*inv(P_eye*(P_eye'))*cur_point1_hat;
    
    L = point_3D_1*point_3D_2' - point_3D_2*point_3D_1';
    cur_point_3D = L * P_i_vec;
    points_3D_homo(:,i) = cur_point_3D;
end