% df/dx, df/du

clear all; clc;
%% dynamics ===========================================================
n_target = 1;
target_dim = 4;
belief_dim = 10;

n_uav = 4;
uav_dim = 4;

% P0 : 센서 파워
% psi : uav 해딩 각도
% beta : 분모 거리의 제곱항
% lambda : 센서와 표적간의 각도 차이
% gamma : 각도 차이 손실 gamma^lambda
% down_anlge : 센서가 아래방향으로 기운 각도
% rotation_angle : 센서가 uav 해딩방향에서 회전한 각도

X = sym('x',[1,(target_dim+belief_dim)*n_target + uav_dim*n_uav],'real');
syms beta coeff_lambda gamma P0 down_angle rotation_angle real
z = sym('z',[1,n_uav],'real'); % 각 uav 고도
down_angle = sym('down_angle',[1,n_uav],'real'); 
rotation_angle = sym('rotation_angle',[1,n_uav],'real');

id_m = 1; 
for id_u = 1:n_uav
    for id_t = 1:n_target
           
        Z_psi2 = [cos(X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+3)), -sin(X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+3)), 0;...
                  sin(X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+3)),  cos(X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+3)), 0;...
                                                                                0,                                                                0, 1];

        X_phi = [1,                                                         0,                                                          0;...
                 0, cos(X(1,(target_dim+belief_dim)*n_target + uav_dim*id_u)), -sin(X(1,(target_dim+belief_dim)*n_target + uav_dim*id_u));...
                 0, sin(X(1,(target_dim+belief_dim)*n_target + uav_dim*id_u)),  cos(X(1,(target_dim+belief_dim)*n_target + uav_dim*id_u))];

        Z_psi1 = [cos(rotation_angle(1,id_u)), -sin(rotation_angle(1,id_u)), 0;...
                  sin(rotation_angle(1,id_u)),  cos(rotation_angle(1,id_u)), 0;...
                                            0,                            0, 1];     

        Y_theta = [ cos(down_angle(1,id_u)), 0,  sin(down_angle(1,id_u));...
                                          0, 1,                        0;...
                   -sin(down_angle(1,id_u)), 0,  cos(down_angle(1,id_u))];               
               
               
        x=[1,0,0]';
        
        unit_vec_sensor = Z_psi2*X_phi*Z_psi1*Y_theta*x; 
        unit_vec_target = [X(1,target_dim*(id_t-1)+1) - X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+1), X(1,target_dim*(id_t-1)+2) - X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+2), -z(1,id_u)]';%/sqrt((X(1,target_dim*(id_t-1)+1) - X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+1))^2+(X(1,target_dim*(id_t-1)+2) - X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1)+2))^2+z(1,id_u)^2);
        unit_vec_target = unit_vec_target/norm(unit_vec_target);
        lambda = acos(unit_vec_target'*unit_vec_sensor);    
        
        sensor_model(id_m,1) = P0*gamma^-(coeff_lambda*lambda^2)/...
        ((X(1,target_dim*(id_t-1) + 1)- X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1) + 1))^2+(X(1,target_dim*(id_t-1) + 2) - X(1,(target_dim+belief_dim)*n_target + uav_dim*(id_u-1) + 2))^2+z(1,id_u)^2)^(beta/2);  
        

        H(id_m,:) = jacobian(sensor_model(id_m,1),X(1,1:target_dim*n_target));
%         jacobian(sensor_model(id_m,1),X(1,1:target_dim*n_target))
        id_m = id_m +1;
    end
end

% for ii = 1:20
%     for jj = 1:40
%         if H(ii,jj) ~= 0
%             ii
%             jj
%             H(ii,jj)
%         end
%     end
% end


% syms dt q R real
% 
% F = [1, 0, dt,  0;...
%      0, 1,  0, dt;...
%      0, 0,  1,  0;...
%      0, 0,  0,  1]; 
% Q =  q*[dt^3/3,      0, dt^2/2,      0;...
%              0, dt^3/3,      0, dt^2/2;...
%         dt^2/2,      0,     dt,      0;...
%              0, dt^2/2,      0,     dt];  
% sigma = [X(1,target_dim*n_target+belief_dim*(id_t-1)+1), X(1,target_dim*n_target+belief_dim*(id_t-1)+2), X(1,target_dim*n_target+belief_dim*(id_t-1)+4), X(1,target_dim*n_target+belief_dim*(id_t-1)+7);...
%          X(1,target_dim*n_target+belief_dim*(id_t-1)+2), X(1,target_dim*n_target+belief_dim*(id_t-1)+3), X(1,target_dim*n_target+belief_dim*(id_t-1)+5), X(1,target_dim*n_target+belief_dim*(id_t-1)+8);...
%          X(1,target_dim*n_target+belief_dim*(id_t-1)+4), X(1,target_dim*n_target+belief_dim*(id_t-1)+5), X(1,target_dim*n_target+belief_dim*(id_t-1)+6), X(1,target_dim*n_target+belief_dim*(id_t-1)+9);...
%          X(1,target_dim*n_target+belief_dim*(id_t-1)+7), X(1,target_dim*n_target+belief_dim*(id_t-1)+8), X(1,target_dim*n_target+belief_dim*(id_t-1)+9), X(1,target_dim*n_target+belief_dim*(id_t-1)+10)];
% 
% x_next = F*X(1,1:target_dim*n_target)';
% 
% sigma_bar = F*sigma*F'+Q;
% K = sigma_bar*H'/(H*sigma_bar*H'+R);
% sigma_new = sigma_bar - K*H*sigma_bar;   
% 
% vec_sigma_new = [sigma_new(1,1), sigma_new(2,1), sigma_new(2,2), sigma_new(3,1), sigma_new(3,2), sigma_new(3,3), sigma_new(4,1), sigma_new(4,2), sigma_new(4,3), sigma_new(4,4)]';
% 
% dfdx = jacobian(vec_sigma_new,X');
% 


















