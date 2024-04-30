%    
% % = sensor model ==========================================================
% 
% clear all; clc;
% syms xt yt xu yu z0 psi beta lambda coeff_lambda gamma P0 down_angle rotation_angle bank_angle real
% 
% % P0 : 센서 파워
% % psi : uav 해딩 각도
% % beta : 분모 거리의 제곱항
% % lambda : 센서와 표적간의 각도 차이
% % gamma : 각도 차이 손실 gamma^lambda
% % down_anlge : 센서가 아래방향으로 기운 각도
% % rotation_angle : 센서가 uav 해딩방향에서 회전한 각도
% % 
% 
% Z_psi2 = [cos(psi), -sin(psi), 0;...
%           sin(psi),  cos(psi), 0;...
%                  0,         0, 1];
% 
% X_phi = [1,               0,                0;...
%          0, cos(bank_angle), -sin(bank_angle);...
%          0, sin(bank_angle),  cos(bank_angle)];
% 
% Z_psi1 = [cos(rotation_angle), -sin(rotation_angle), 0;...
%           sin(rotation_angle),  cos(rotation_angle), 0;...
%                                0,                        0, 1];     
%  
% Y_theta = [ cos(down_angle), 0,  sin(down_angle);...
%                           0, 1,                0;...
%            -sin(down_angle), 0,  cos(down_angle)];
%        
%        
%        
% x=[1,0,0]';
% unit_vec_sensor = Z_psi2*X_phi*Z_psi1*Y_theta*x;   
% 
% unit_vec_target = [xt-xu, yt-yu, -z0]'/sqrt((xt-xu)^2+(yt-yu)^2+z0^2);
% 
% lambda = acos(unit_vec_target'*unit_vec_sensor);
% 
% z = P0*gamma^-(coeff_lambda*lambda^2)/((xt-xu)^2+(yt-yu)^2+z0^2)^(beta/2);   
% 
% % 
% % H11 = diff(z,xt)
% % H21 = diff(z,yt)


% = final test ============================================================
clear all; clc; 

alpha = 1000;
z0 = 5; % m
beta = 100;
coeff_lambda = 3;
gamma = 1.5;

down_angle = (30)/180*pi;     % y-axis   아래가 +  
bank_angle = (0)/180*pi;
        
rotation_angle = (0)/180*pi; % z-axis  회전 방향(uav 해딩 방향+)
psi = (0)/180*pi;

Z_psi2 = [cos(psi), -sin(psi), 0;...
          sin(psi),  cos(psi), 0;...
                 0,         0, 1];

X_phi = [1,               0,                0;...
         0, cos(bank_angle), -sin(bank_angle);...
         0, sin(bank_angle),  cos(bank_angle)];

Z_psi1 = [cos(rotation_angle), -sin(rotation_angle), 0;...
          sin(rotation_angle),  cos(rotation_angle), 0;...
                            0,                    0, 1];     
 
Y_theta = [ cos(down_angle), 0,  sin(down_angle);...
                          0, 1,                0;...
           -sin(down_angle), 0,  cos(down_angle)];
     
x=[1,0,0]';
% unit_vec_sensor = Z_psi2*X_phi*Z_psi1*Y_theta*x; 
unit_vec_sensor = [1,0]';

xu = 0;
yu = 0;
xxt = -50:0.2:50;
yyt = -50:0.2:50;

[X, Y] = meshgrid(xxt, yyt);

figure; hold on; grid on; %zlim([0 200]);         
for i = 1:length(xxt)
    xt = xxt(i);
    for j = 1:length(yyt)
        yt = yyt(j);
        
%         unit_vec_target = [xt-xu, yt-yu, -z0]'/sqrt((xt-xu)^2+(yt-yu)^2+z0^2);
        unit_vec_target = [xt-xu, yt-yu]'/sqrt((xt-xu)^2+(yt-yu)^2);
        lambda = acos(unit_vec_target'*unit_vec_sensor);
        
        z = alpha*gamma^(-coeff_lambda*lambda^2)/(((xt-xu)^2+(yt-yu)^2+z0^2)^(2/2)+beta);   
%         z = P0/(gamma^(coeff_lambda*acos(((yt - yu)*(cos(down_angle)*(cos(rotation_angle)*sin(psi) + cos(bank_angle)*cos(psi)*sin(rotation_angle)) + cos(psi)*sin(bank_angle)*sin(down_angle)))/((xt - xu)^2 + (yt - yu)^2 + z0^2)^(1/2) + ((xt - xu)*(cos(down_angle)*(cos(psi)*cos(rotation_angle) - cos(bank_angle)*sin(psi)*sin(rotation_angle)) - sin(bank_angle)*sin(down_angle)*sin(psi)))/((xt - xu)^2 + (yt - yu)^2 + z0^2)^(1/2) + (z0*(cos(bank_angle)*sin(down_angle) - cos(down_angle)*sin(bank_angle)*sin(rotation_angle)))/((xt - xu)^2 + (yt - yu)^2 + z0^2)^(1/2))^2)*((xt - xu)^2 + (yt - yu)^2 + z0^2)^(beta/2));

        Z(j,i) = z;
    end    
end
mesh(X,Y,Z)
plot3(xu,yu,1,'ro')
xlabel('X')
ylabel('Y')
view([0 90])

max(max(Z))
z0/tan(down_angle)
