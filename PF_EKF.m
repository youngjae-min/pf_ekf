clear; clc;

nv = 4;         % number of agents (vehicles)
np = 500;       % number of particles
dt = 1;         % time interval for each step
N = 50;         % number of points to use in Gauss Quadrature
maxiter = 1;    % maximum #iteration for alternating optimization
hoffmann = 1;   % 1: run hoffmann 0: don't

%% motion model
% fixed-wing UAV model
xdim = 3;       % dimension of agent's state [x-coord y-coord heading-ang bank-ang]
udim = 2;       % dimension of control input (bank_change)
v = 1;          % constant speed [m/s]
ga = 9.8;       % gravitational acceleration [m/s^2]
Rmin = 3;     % minimum radius of gyration that UAV can acheive [m] 
P = [0.01    0 0;...   % covariance matrix for noise
     0    0.01 0;...
     0      0   0.0076];
P = 0.5*P;
% P = zeros(3);
num_action = 8; % number of actions to consider
Phimax = atan(v^2/Rmin/ga); % limit on roll angle
action_set1 = v;%linspace(0.01,2,num_action);
action_set2 = linspace(-Phimax,Phimax,num_action);

p_s = sym('p',[xdim 1], 'real');
u_s = sym('u',[udim 1],'real');
h_s = p_s + [u_s(1)*cos(p_s(3))*dt; u_s(1)*sin(p_s(3))*dt; ga*tan(u_s(2))*dt/u_s(1)];

% % directional movement
% xdim = 2;       % dimension of agent's state [x-coord y-coord]
% udim = 1;       % dimension of control input (direction)
% v = 1;          % constant speed [m/s]
% num_action = 8; % number of actions to consider
% P = [0.6   0.4 ;...   % covariance matrix for noise
%      0.4   0.6]; %zeros(2);
% action_set = (0:num_action-1)*2*pi/num_action-pi;         % limit control input in from 0 to 2pi
% 
% p_s = sym('p',[xdim 1],'real');
% syms u_s real
% h_s = p_s + [v*cos(u_s)*dt; v*sin(u_s)*dt];

% % double integrator
% xdim = 4;       % dimension of agent's state [x-coord y-coord x-speed y-speed]
% udim = 2;       % dimension of control input [x-accel y-accel]
% P = [0.01    0 0   0;...   % covariance matrix for noise
%      0      0.01 0   0;...
%      0      0   0.0001 0;...
%      0      0   0   0.0001];
% % P = zeros(4);
% vlim = 1;       % speed limit [m/s]
% num_action = 8; % number of actions to consider
% ulim = 0.2;         % limit for control input (acceleration)
% action_set1 = linspace(-ulim,ulim,num_action);
% action_set2 = linspace(-ulim,ulim,num_action);
% 
% p_s = sym('p',[xdim 1],'real');
% u_s = sym('u',[udim 1],'real');
% h_s = p_s + [p_s(3)*dt+u_s(1)*dt^2/2; p_s(4)*dt+u_s(2)*dt^2/2; u_s(1)*dt; u_s(2)*dt];

% % vel-ang model
% xdim = 3;       % dimension of agent's state [x-coord y-coord angle]
% udim = 2;       % dimension of control input [speed ang_vel]
% P = [0.1    0   0;...   % covariance matrix for noise
%      0      0.1 0;...
%      0      0   0.01];
% % P = zeros(3);
% vlim = 1;       % speed limit [m/s]
% num_action = 8; % number of actions to consider
% action_set1 = linspace(0,2,num_action);
% action_set2 = linspace(-pi/6,pi/6,num_action);
% p_s = sym('p',[xdim 1],'real');
% u_s = sym('u',[udim 1],'real');
% h_s = p_s + [u_s(1)/u_s(2)*(sin(p_s(3)+u_s(2)*dt)-sin(p_s(3))); u_s(1)/u_s(2)*(cos(p_s(3))-cos(p_s(3)+u_s(2)*dt)); u_s(2)*dt];

% vpa(h_s,7);       % integer expression to floating-point arithmetic
H_s = jacobian(h_s,p_s);

% change to function
h_func = matlabFunction(h_s,'Vars',{p_s,u_s});
H_func = matlabFunction(H_s,'Vars',{p_s,u_s});

%% observation model
pt_s = sym('pt', [2 1], 'real');  % target's state [x-coord y-coord]
z0 = 5;                 % fixed height of agents [m]
R = 2;                  % noise covariance
% zmax = 121.4939;% 200;        % maximum observation value
% zmin = 0;               % minimum observation value

% SNR
alpha = 1000;       
beta = 100;
coeff_lambda = 3;
gamma = 1.5;

% % sensor direction reference
% down_angle = (30)/180*pi;       % y-axis   아래가 +         
% rotation_angle = (0)/180*pi;    % z-axis  회전 방향(uav 해딩 방향+)
% 
% Z_psi2 = [cos(p_s(3)), -sin(p_s(3)), 0;...
%           sin(p_s(3)),  cos(p_s(3)), 0;...
%                   0,          0, 1];
% 
% X_phi = [1,         0,          0;...
%          0, cos(p_s(4)), -sin(p_s(4));...
%          0, sin(p_s(4)),  cos(p_s(4))];
% 
% Z_psi1 = [cos(rotation_angle), -sin(rotation_angle), 0;...
%           sin(rotation_angle),  cos(rotation_angle), 0;...
%                             0,                    0, 1];     
%  
% Y_theta = [ cos(down_angle), 0,  sin(down_angle);...
%                           0, 1,                0;...
%            -sin(down_angle), 0,  cos(down_angle)];
%      
% ref = [1; 0; 0];
% unit_vec_sensor = Z_psi2*X_phi*Z_psi1*Y_theta*ref; 
% unit_vec_target = [p_s(1)-pt_s(1); p_s(2)-pt_s(2); -z0];
% unit_vec_target = unit_vec_target/norm(unit_vec_target);
% lambda = acos(unit_vec_target'*unit_vec_sensor);    % angle difference between sensor and target
% 
% g = P0*gamma^-(coeff_lambda*lambda^2)/((p_s(1)-pt_s(1))^2+(p_s(2)-pt_s(2))^2+z0^2)^(beta/2);

unit_vec_sensor = [cos(p_s(3)); sin(p_s(3))]; 
unit_vec_target = [pt_s(1)-p_s(1); pt_s(2)-p_s(2)];
unit_vec_target = unit_vec_target/norm(unit_vec_target);
lambda = acos(unit_vec_target'*unit_vec_sensor);
g_s = alpha*gamma^-(coeff_lambda*lambda^2)/(((p_s(1)-pt_s(1))^2+(p_s(2)-pt_s(2))^2+z0^2)^(2/2)+beta);
% g_s = lambda;

% % distnace
% g_s = 2000/(sqrt((p_s(1)-pt_s(1))^2+(p_s(2)-pt_s(2))^2+z0^2)^2+100);

% vpa(g_s,7);       % integer expression to floating-point arithmetic
G_s = jacobian(g_s,p_s);

% change to function
g_func = matlabFunction(g_s,'Vars',{p_s,pt_s});
G_func = matlabFunction(G_s,'Vars',{p_s,pt_s});

%% visualization setting
scaling = 2.4477;           % scaling factor for 95% confidence ellipse
t = linspace(0,2*pi,20);    % 20 points for drawing ellipse
e = [cos(t) ; sin(t)];      % unit circle
color = {[0.4940 0.1840 0.5560],[0.8500 0.3250 0.0980],[0.3010 0.7450 0.9330],'m'};  % color for each agent

%% initial position in map [0,40]x[0,40]
% theta_true = [15; 23];
theta_true = [10; 30];%40*rand(2,1);
% x = [[37; 21.5; pi; 0] [37; 18.5; pi; 0] [40; 21.5; pi; 0] [40; 18.5; pi; 0]];
% x = [[37; 21.5] [37; 18.5] [40; 21.5] [40; 18.5]];
% x = [[37; 21.5; 0; 0] [37; 18.5; 0; 0] [40; 21.5; 0; 0] [40; 18.5; 0; 0]];
% x = [[37; 21.5; pi] [37; 18.5; pi] [40; 21.5; pi] [40; 18.5; pi]];
x = [[57; 21.5; pi] [57; 18.5; pi] [60; 21.5; pi] [60; 18.5; pi]];

%% initialize filters
w = ones(1,np)/np;                          % weights of PF
theta = rand(2,np)*60;                      % particles
mu = repmat(reshape(x,xdim,nv,1),1,1,np);   % mu_k^(j) is mu(:,j,k)
sigma = zeros(xdim,xdim,nv,np);             % sigma_k^(j) is sigma(:,:,j,k)

%% main
f = figure();
x = reshape(x,xdim,nv,1);       % record every position of each step along the 3rd axis
x_naive = x;                    % when noise in motion model is not considered
tstamp = 1;                     % time stamp (1 is for initial)
u_max = zeros(udim,nv)+0.001;         % optimized control input
z = zeros(nv,1);                % observation at each step

% % gaussian approx _ gradient descent
% mu_s = sym('mu',[xdim,nv,np],'real');
% sigma_s = sym('sigma',[xdim,xdim,nv,np],'real');
% u_all_s = sym('u_all',[4,1],'real');
% mu_bar_s = sym('mu_bar',[xdim,nv,np],'real');
% sigma_bar_s = sym('sigma_bar',[xdim,xdim,nv,np],'real');
% mu_tilda_s =  sym('mu_tilda',[nv,np],'real');
% sigma_tilda_s = sym('sigma_tilda',[nv,np],'real');
% G_all_s = sym('G_all',[1,xdim,nv,np],'real');
% mu_total_s = sym('mu_total',[nv,1],'real');
% sigma_total_s = sym('sigma_total',[nv,nv],'real');
% for j=1:nv
%     for k=1:np
%         H = H_func(mu_s(:,j,k), u_all_s(j));
%         mu_bar_s(:,j,k) = h_func(mu_s(:,j,k), u_all_s(j));
%         sigma_bar_s(:,:,j,k) = H*sigma_s(:,:,j,k)*H'+P;
%         G_all_s(1,:,j,k) = G_func(mu_s(:,j,k), theta(:,k));
%         mu_tilda_s(j,k) = g_func(mu_bar_s(:,j,k), theta(:,k));
%         sigma_tilda_s(j,k) = G_all_s(1,:,j,k)*sigma_bar_s(:,:,j,k)*G_all_s(1,:,k)'+R;
%     end
%     mu_total_s(j) = w*mu_tilda_s(j,:)';
%     sigma_total_s(j,j) = w*(sigma_tilda_s(j,:)+mu_tilda_s(j,:).^2)' - mu_total_s(j)^2;
%     for i=1:j-1
%         sigma_total_s(i,j) = w*((mu_tilda_s(i,:)-mu_total_s(i)).*((mu_tilda_s(j,:)-mu_total_s(j))))';
%     end
%     for i=j+1:nv
%         sigma_total_s(i,j) = sigma_total_s(j,i);
%     end
% end
% % I_s = log(det(sigma_total_s))- sum(w*log(sigma_tilda_s'));
% I_s = sum(log(diag(GaussElimination(sigma_total_s,''))))- sum(w*log(sigma_tilda_s'));
% grad_s = jacobian(I_s,u_all_s);
% 
% % to matlab function
% grad_func = matlabFunction(grad_s,'Vars',{mu_s,sigma_s,u_all_s});
% mu_bar_func = matlabFunction(mu_bar_s,'Vars',{mu_s,sigma_s,u_all_s});
% sigma_bar_func = matlabFunction(sigma_bar_s,'Vars',{mu_s,sigma_s,u_all_s});
% mu_tilda_func = matlabFunction(mu_tilda_s,'Vars',{mu_s,sigma_s,u_all_s});
% sigma_tilda_func = matlabFunction(sigma_tilda_s,'Vars',{mu_s,sigma_s,u_all_s});
% G_all_func = matlabFunction(G_all_s,'Vars',{mu_s,sigma_s,u_all_s});
% 
% save('functions.mat',grad_func,mu_bar_func,sigma_bar_func,mu_tilda_func,sigma_tilda_func,G_all_func);

% gaussian approx _ grid search
mu_bar = zeros(xdim,nv,np);
sigma_bar = zeros(xdim,xdim,nv,np);
mu_tilda =  zeros(nv,np);
sigma_tilda = zeros(nv,np);
G = zeros(1,xdim,nv,np);
mu_total = zeros(nv,1);
sigma_total = zeros(nv,nv);

% % single node approx
% mu_bar = zeros(xdim,np);            % reuse for all j=1:nv
% sigma_bar = zeros(xdim,xdim,np);    % reuse
% mu_tilda =  zeros(1,np);              % reuse
% sigma_tilda = zeros(1,np);            % reuse
% G = zeros(1,xdim,np);                 % reuse
% mu_bar_max = zeros(xdim,nv,np);
% sigma_bar_max = zeros(xdim,xdim,nv,np);
% mu_tilda_max = zeros(nv,np);
% sigma_tilda_max = zeros(nv,np);
% G_max = zeros(1,xdim,nv,np);
% % [zq, wq] = lgwt(N,zmin, zmax);  % Gauss Quadrature points and weights
% zq_raw = randn(N,1);

if hoffmann
    set(gcf, 'Position', [500 300 1000 400]);
    x_h = x;
    x_naive_h = x;
    w_h = w;
    theta_h = theta;
    mu_h = mu_tilda;
    u_max_h = u_max;
    z_h = z;
    mu_total_h = mu_total;
    sigma_total_h = sigma_total;
end

tic
while true
    tstamp = tstamp + 1;
    %% visualize
    if mod(tstamp,10)==0
        clf(f)
        if hoffmann
            subplot(1,2,2)
            axis([-10 70 -10 70])
            hold on
            % PF
            w_max = max(w_h);
            for k=1:np
                scatter(theta_h(1,k), theta_h(2,k),18,'MarkerFaceColor','b',...
                    'MarkerEdgeColor','none','MarkerFaceAlpha', 0.05+0.95*w_h(k)/w_max)
            end
            % target and agents
            theta_mmse_h = theta_h*w_h';
            scatter(theta_true(1), theta_true(2), 50, 'filled', 'sr')
            scatter(theta_mmse_h(1), theta_mmse_h(2), 50, 'xg','LineWidth',2)
            for j=1:nv
                plot(squeeze(x_h(1,j,:)), squeeze(x_h(2,j,:)),'LineStyle',':','Color',color{j},'LineWidth',2)
                scatter(x_h(1,j,end), x_h(2,j,end),200,'p','LineWidth',1,'MarkerFaceColor',color{j},'MarkerEdgeColor','w')
                scatter(x_naive_h(1,j,end), x_naive_h(2,j,end),150,'p','LineWidth',0.8,'MarkerFaceColor','none','MarkerEdgeColor',color{j})
            end
            hold off
            subplot(1,2,1)
        end
        axis([-10 70 -10 70])
        hold on
        % PF
        w_max = max(w);
        for k=1:np
            scatter(theta(1,k), theta(2,k),18,'MarkerFaceColor','b',...
                'MarkerEdgeColor','none','MarkerFaceAlpha', 0.05+0.95*w(k)/w_max)
        end
    %     % EKF - separated
    %     for j=1:nv
    %         for k=1:np
    %             if w(k)>w_max/2
    %                 % eigen decomposition [sorted by eigen values]
    %                 [V, D] = eig(sigma(1:2,1:2,j,k)); % consider position only
    %                 [D, order] = sort(diag(D), 'descend');
    %                 D = diag(D);
    %                 V = V(:, order); 
    %                 ell = 2.447*V*sqrt(D)*e + mu(1:2,j,k);  % scale eigenvectors and project circle back to orig space
    %                 plot(ell(1,:), ell(2,:), 'Color',color{j}); % plot ellipse
    %             end
    %         end
    %     end
        % target and agents
        theta_mmse = theta*w';
        scatter(theta_true(1), theta_true(2), 50, 'filled', 'sr')
        scatter(theta_mmse(1), theta_mmse(2), 50, 'xg','LineWidth',2)
        for j=1:nv
            plot(squeeze(x(1,j,:)), squeeze(x(2,j,:)),'LineStyle',':','Color',color{j},'LineWidth',2)
            scatter(x(1,j,end), x(2,j,end),200,'p','LineWidth',1,'MarkerFaceColor',color{j},'MarkerEdgeColor','w')
    %         scatter(x_naive(1,j,end), x_naive(2,j,end),150,'p','LineWidth',0.8,'MarkerFaceColor','none','MarkerEdgeColor',color{j})
        end
         % EKF - integrated
        for j=1:nv
            mu_int = squeeze(mu(1:2,j,:))*w';
            sigma_int = -mu_int*mu_int';
            for k=1:np
                sigma_int = sigma_int + w(k)*(sigma(1:2,1:2,j,k)+mu(1:2,j,k)*mu(1:2,j,k)');
            end
            [V, D] = eig(sigma_int); % consider position only
            [D, order] = sort(diag(D), 'descend');
            D = diag(D);
            V = V(:, order); 
            ell = 2.447*V*sqrt(D)*e + mu_int;  % scale eigenvectors and project circle back to orig space
            plot(ell(1,:), ell(2,:), 'Color',color{j}); % plot ellipse
        end
        hold off
        drawnow
        disp('visualize done')
        toc
    end
    %% planning
%     % gaussian approx _ grid search
%     I_max = -inf;
%     for u1=action_set
%         for u2=action_set
%             for u3=action_set
%                 for u4=action_set
%                     u = [u1 u2 u3 u4];
%                     for j=1:nv
%                         for k=1:np
%                             H = H_func(mu(:,j,k), u(j));
%                             mu_bar(:,j,k) = h_func(mu(:,j,k), u(j));
%                             sigma_bar(:,:,j,k) = H*sigma(:,:,j,k)*H'+P;
%                             G(1,:,j,k) = G_func(mu(:,j,k), theta(:,k));
%                             mu_tilda(j,k) = g_func(mu_bar(:,j,k), theta(:,k));
%                             sigma_tilda(j,k) = G(1,:,j,k)*sigma_bar(:,:,j,k)*G(1,:,j,k)'+R;
%                         end
%                         mu_total(j) = w*mu_tilda(j,:)';
%                         sigma_total(j,j) = w*(sigma_tilda(j,:)+mu_tilda(j,:).^2)' - mu_total(j)^2;
%                         for i=1:j-1
%                             sigma_total(i,j) = w*((mu_tilda(i,:)-mu_total(i)).*((mu_tilda(j,:)-mu_total(j))))';
%                             sigma_total(j,i) = sigma_total(i,j);
%                         end
%                     end
%                     
%                     I = log(det(sigma_total))- sum(w*log(sigma_tilda'));
%                     
%                     if I>I_max
%                         I_max = I;
%                         u_max = u;
%                         mu_bar_max = mu_bar;
%                         sigma_bar_max = sigma_bar;
%                         mu_tilda_max = mu_tilda;
%                         sigma_tilda_max = sigma_tilda;
%                         G_max = G;
%                     end
%                 end
%             end
%         end
%     end

    % gaussian approx _ alternating grid search / double integrator
    for iter=1:maxiter
        for jj=1:nv
            I_max = -inf;
            u = u_max;
            for u1=action_set1
                for u2=action_set2
                    u(:,jj) = [u1; u2];
                    if iter==1 && jj==1 && I_max == -inf % set all things to new
                        for j=1:nv
                            for k=1:np
                                H = H_func(mu(:,j,k), u(:,j));
                                mu_bar(:,j,k) = h_func(mu(:,j,k), u(:,j));
                                sigma_bar(:,:,j,k) = H*sigma(:,:,j,k)*H'+P;
                                G(1,:,j,k) = G_func(mu(:,j,k), theta(:,k));
                                mu_tilda(j,k) = g_func(mu_bar(:,j,k), theta(:,k));
                                sigma_tilda(j,k) = G(1,:,j,k)*sigma_bar(:,:,j,k)*G(1,:,j,k)'+R;
                            end
                            mu_total(j) = w*mu_tilda(j,:)';
                            sigma_total(j,j) = w*(sigma_tilda(j,:)+mu_tilda(j,:).^2)' - mu_total(j)^2;
                            for i=1:j-1
                                sigma_total(i,j) = w*((mu_tilda(i,:)-mu_total(i)).*((mu_tilda(j,:)-mu_total(j))))';
                                sigma_total(j,i) = sigma_total(i,j);
                            end
                        end
                    else          % change only what changes
                        for k=1:np
                            H = H_func(mu(:,jj,k), u(:,jj));
                            mu_bar(:,jj,k) = h_func(mu(:,jj,k), u(:,jj));
%                             if max(abs(mu_bar(:,jj,k))) > 150
%                                 disp('stop')
%                             end
                                
                            sigma_bar(:,:,jj,k) = H*sigma(:,:,jj,k)*H'+P;
                            G(1,:,jj,k) = G_func(mu(:,jj,k), theta(:,k));
                            mu_tilda(jj,k) = g_func(mu_bar(:,jj,k), theta(:,k));
                            sigma_tilda(jj,k) = G(1,:,jj,k)*sigma_bar(:,:,jj,k)*G(1,:,jj,k)'+R;
                        end
                        mu_total_old = mu_total(jj);
                        sigma_total_old = sigma_total(jj,:);
                        mu_total(jj) = w*mu_tilda(jj,:)';
                        sigma_total(jj,jj) = w*(sigma_tilda(jj,:)+mu_tilda(jj,:).^2)' - mu_total(jj)^2;
                        for i=[1:jj-1 jj+1:nv]
                            sigma_total(i,jj) = w*((mu_tilda(i,:)-mu_total(i)).*((mu_tilda(jj,:)-mu_total(jj))))';
                            sigma_total(jj,i) = sigma_total(i,jj);
                        end
                    end

                    I = log(det(sigma_total))- sum(w*log(sigma_tilda'));

%                     if max(abs(mu_bar(:,jj,1))) > 150
%                         disp('stop')
%                     end
                            
                    if I>I_max
                        I_max = I;
                        u_max = u;
                        mu_bar_max = mu_bar;
                        sigma_bar_max = sigma_bar;
                        mu_tilda_max = mu_tilda;
                        sigma_tilda_max = sigma_tilda;
                        G_max = G;
                    else    % undo the invalid changes
                        mu_bar(:,jj,:) = mu_bar_max(:,jj,:);
                        sigma_bar(:,:,jj,:) = sigma_bar_max(:,:,jj,:);
                        G(1,:,jj,:) = G_max(1,:,jj,:);
                        mu_tilda(jj,:) = mu_tilda_max(jj,:);
                        sigma_tilda(jj,:) = sigma_tilda_max(jj,:);
                        mu_total(jj) = mu_total_old;
                        sigma_total(jj,:) = sigma_total_old;
                        sigma_total(:,jj) = sigma_total(jj,:)';
                    end
                end
            end
        end
    end
    
    if hoffmann
%         u_max_h = u_max;
%         u_h = u_max;
%         for j=1:nv
%             for k=1:np
%                 mu_h(j,k) = g_func(h_func(x_h(:,j,tstamp-1), u_h(:,j)), theta_h(:,k));
%             end
%         end
        for iter=1:maxiter
            for jj=1:nv
                I_max = -inf;
                u_h = u_max_h;
                for u1=action_set1
                    for u2=action_set2
                        u_h(:,jj) = [u1; u2];
                        if iter==1 && jj==1 && I_max == -inf % set all things to new
                            for j=1:nv
                                for k=1:np
                                    mu_h(j,k) = g_func(h_func(x_h(:,j,tstamp-1), u_h(:,j)), theta_h(:,k));
                                end
                                mu_total_h(j) = w_h*mu_h(j,:)';
                                sigma_total_h(j,j) = R+w_h*(mu_h(j,:).^2)' - mu_total_h(j)^2;
                                for i=1:j-1
                                    sigma_total_h(i,j) = w_h*((mu_h(i,:)-mu_total_h(i)).*((mu_h(j,:)-mu_total_h(j))))';
                                    sigma_total_h(j,i) = sigma_total_h(i,j);
                                end
                            end
                        else          % change only what changes
                            mu_old = mu_h(jj,:);
                            mu_total_old = mu_total_h(jj);
                            sigma_total_old = sigma_total_h(jj,:);
                            for k=1:np
                                mu_h(jj,k) = g_func(h_func(x_naive_h(:,jj,tstamp-1), u_h(:,jj)), theta_h(:,k));
                            end
                            mu_total_h(jj) = w_h*mu_h(jj,:)';
                            sigma_total_h(jj,jj) = R+w_h*(mu_h(jj,:).^2)' - mu_total_h(jj)^2;
                            for i=[1:jj-1 jj+1:nv]
                                sigma_total_h(i,jj) = w_h*((mu_h(i,:)-mu_total_h(i)).*((mu_h(jj,:)-mu_total_h(jj))))';
                                sigma_total_h(jj,i) = sigma_total_h(i,jj);
                            end
                        end
                        
                        I = log(det(sigma_total_h));

                        if I>I_max
                            I_max = I;
                            u_max_h = u_h;
                        else
                            mu_h(jj,:) = mu_old;
                            mu_total_h(jj) = mu_total_old;
                            sigma_total_h(jj,:) = sigma_total_old;
                            sigma_total_h(:,jj) = sigma_total_h(jj,:)';
                        end
                    end
                end
            end
        end
    end
    
%     % single node approx / double integrator
%     for j=1:nv
%         I_max = -inf;
%         v_avg = squeeze(mu(3:4,j,:))*w';
%         for u1=action_set
%             for u2=action_set
%                 u=[u1; u2];
%                 if u'*u > ulim^2  % over the acceleration limit
%                     continue
%                 end
%                 if sum((v_avg+u*dt).^2) > vlim^2   % over the speed limit
%                     continue
%                 end
%                 for k=1:np
%                     H = H_func(mu(:,j,k), u);
%                     mu_bar(:,k) = h_func(mu(:,j,k), u);
%                     sigma_bar(:,:,k) = H*sigma(:,:,j,k)*H'+P;
%                     G(1,:,k) = G_func(mu(:,j,k), theta(:,k));
%                     mu_tilda(k) = g_func(mu_bar(:,k), theta(:,k));
%                     sigma_tilda(k) = G(1,:,k)*sigma_bar(:,:,k)*G(1,:,k)'+R;
%                 end
%                 sig = repmat(sigma_tilda,N,1);
%     %             I = exp(-(repmat(zq,1,np)-repmat(mu_tilda,N,1)).^2./sig/2)./sqrt(sig);
%     %             I = I*w';
%     %             I = -wq'*(I.*log(I)) - w*log(sigma_tilda)'/2;
%                 I = - w*log(sigma_tilda)'/2;
%                 for k=1:np
%                     zq = mu_tilda(k) + sqrt(sigma_tilda(k))*zq_raw;
%                     Ik = exp(-(repmat(zq,1,np)-repmat(mu_tilda,N,1)).^2./sig/2)./sqrt(sig);
%                     I = I- w(k)*sum(log(nonzeros(Ik*w')));
%     %                 if isnan(I)
%     %                     disp('here!')
%     %                 end
%                 end
% 
%                 if I>I_max
%                     I_max = I;
%                     u_max(:,j) = u;
%                     mu_bar_max(:,j,:) = mu_bar;
%                     sigma_bar_max(:,:,j,:) = sigma_bar;
%                     mu_tilda_max(j,:) = mu_tilda;
%                     sigma_tilda_max(j,:) = sigma_tilda;
%                     G_max(1,:,j,:) = G;
%                 end
%             end
%         end
%     end

%     % single node approx
%     for j=1:nv
%         I_max = -inf;
% %         bank_avg = w*squeeze(mu(4,j,:));
%         for u=action_set
% %             if abs(bank_avg+u_*dt) > Phimax % beyond the limit
% %                 continue
% %             end
%             for k=1:np
%                 H = H_func(mu(:,j,k), u);
%                 mu_bar(:,k) = h_func(mu(:,j,k), u);
%                 sigma_bar(:,:,k) = H*sigma(:,:,j,k)*H'+P;
%                 G(1,:,k) = G_func(mu(:,j,k), theta(:,k));
%                 mu_tilda(k) = g_func(mu_bar(:,k), theta(:,k));
%                 sigma_tilda(k) = G(1,:,k)*sigma_bar(:,:,k)*G(1,:,k)'+R;
%             end
%             sig = repmat(sigma_tilda,N,1);
% %             I = exp(-(repmat(zq,1,np)-repmat(mu_tilda,N,1)).^2./sig/2)./sqrt(sig);
% %             I = I*w';
% %             I = -wq'*(I.*log(I)) - w*log(sigma_tilda)'/2;
%             I = - w*log(sigma_tilda)'/2;
%             for k=1:np
%                 zq = mu_tilda(k) + sqrt(sigma_tilda(k))*zq_raw;
%                 Ik = exp(-(repmat(zq,1,np)-repmat(mu_tilda,N,1)).^2./sig/2)./sqrt(sig);
%                 I = I- w(k)*sum(log(nonzeros(Ik*w')));
% %                 if isnan(I)
% %                     disp('here!')
% %                 end
%             end
% 
%             if I>I_max
%                 I_max = I;
%                 u_max(j) = u;
%                 mu_bar_max(:,j,:) = mu_bar;
%                 sigma_bar_max(:,:,j,:) = sigma_bar;
%                 mu_tilda_max(j,:) = mu_tilda;
%                 sigma_tilda_max(j,:) = sigma_tilda;
%                 G_max(1,:,j,:) = G;
%             end
%         end
%     end
    disp('planning done')
    toc
    %% move and observe
    for j=1:nv
        noise_motion = mvnrnd(zeros(xdim,1), P)';
        noise_observ = randn(1)*sqrt(R);
        if hoffmann
            x_h(:,j,tstamp) = h_func(x_h(:,j,tstamp-1), u_max_h(:,j)) + noise_motion;
            x_naive_h(:,j,tstamp) = h_func(x_naive_h(:,j,tstamp-1), u_max_h(:,j));
            z_h(j) = g_func(x_h(:,j,tstamp-1), theta_true) + noise_observ;
        end
        x(:,j,tstamp) = h_func(x(:,j,tstamp-1), u_max(:,j)) + noise_motion;
%         x_naive(:,j,tstamp) = h_func(x_naive(:,j,tstamp-1), u_max(:,j));
        z(j) = g_func(x(:,j,tstamp-1), theta_true) + noise_observ;
    end
    disp('move&observe done')
    toc
    %% update
    r = rand(1)/np; % random seed
    
    for k=1:np
        for j=1:nv
            K = sigma_bar_max(:,:,j,k) * G_max(1,:,j,k)' /sigma_tilda_max(j,k);
            mu(:,j,k) = mu_bar_max(:,j,k) + K*(z(j)-mu_tilda_max(j,k));
%             if max(abs(mu(:,j,k))) > 150
%                 disp('stop')
%             end
            sigma(:,:,j,k) = sigma_bar_max(:,:,j,k) - K*G_max(1,:,j,k)*sigma_bar_max(:,:,j,k);
            w(k) = w(k)*exp(-(z(j)-mu_tilda_max(j,k))^2/(2*sigma_tilda_max(j,k)))/sqrt(sigma_tilda_max(j,k));
        end
    end
    w = w/sum(w);   % normalize
    % PF resampling
    if w*w'>2/np % when Neff < np/2
        disp('doing resampling..')
        theta_old = theta;
        mu_old = mu;
        sigma_old = sigma;
%         r = rand(1)/np;
        c = w(1);
        i = 1;
        for m=1:np
            U = r+(m-1)/np;
            while U>c
                i = i+1;
                c = c+w(i);
            end
            theta(:,m) = theta_old(:,i);
            mu(:,:,m) = mu_old(:,:,i);
            sigma(:,:,:,m) = sigma_old(:,:,:,i);
        end
        theta = theta + 0.3*randn(2,np);
        w = ones(1,np)/np;
    end
    
    if hoffmann
        w_h = w_h.*prod(exp(-(repmat(z_h,1,np)-mu_h).^2/(2*R)),1);
        w_h = w_h/sum(w_h);   % normalize
        
        % PF resampling
        if w_h*w_h'>2/np % when Neff < np/2
            disp('doing resampling.. (hoffmann)')
            theta_old = theta_h;
            c = w_h(1);
            i = 1;
            for m=1:np
                U = r+(m-1)/np;
                while U>c
                    i = i+1;
                    c = c+w_h(i);
                end
                theta_h(:,m) = theta_old(:,i);
            end
            theta_h = theta_h + 0.3*randn(2,np);
            w_h = ones(1,np)/np;
        end
    end
    
    disp('done update')
    toc
%     disp('click for next step')
%     waitforbuttonpress;
end     