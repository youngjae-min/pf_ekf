clear all;
close all;

% eigen decomposition [sorted by eigen values]
[V, D] = eig(covariance);
[D, order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);

t = linspace(0,2*pi,20);
e = [cos(t) ; sin(t)];        % unit circle
e = 2.447*V*sqrt(D)*e + Mu';  % scale eigenvectors and project circle back to orig space

% plot cov and major/minor axes
plot(e(1,:), e(2,:), 'Color','k');
