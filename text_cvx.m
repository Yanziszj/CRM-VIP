% Clear workspace
clear; clc;

% Problem size
n = 2;   % Dimension
m = 1;   % Number of ellipsoids

% Point to project
x0 = randn(n, 1);
x0 = [100;0];

% Store ellipsoid parameters
A = cell(m, 1);  % symmetric matrices
b = cell(m, 1);  % vectors
c = zeros(m, 1); % scalars

% Random ellipsoids (make sure A{i} is symmetric positive semidefinite)
for i = 1:m
    Q = randn(n, n);
    A{i} = Q' * Q + 1e-3 * eye(n);  % ensure PSD
    b{i} = randn(n, 1);
    c(i) = rand() * 5;              % some positive scalar
end

% Use CVX to compute projection
cvx_begin quiet
    variable x(n)
    minimize( sum_square(x - x0) )
    subject to
        for i = 1:m
            quad_form(x, A{i}) + 2 * b{i}' * x - c(i) <= 0;
        end
cvx_end

% Output
disp('Original point x0:');
disp(x0);
disp('Projected point x:');
disp(x);

% Check constraint values
fprintf('\nConstraint check (should be ¡Ü 0):\n');
for i = 1:m
    val = x' * A{i} * x + 2 * b{i}' * x - c(i);
    fprintf('Ellipsoid %d: %.6f\n', i, val);
end

 (x-x0)'*(A{1}*x+2*b{i})/(norm(x-x0)*norm(A{1}*x+2*b{i}))


C=cell(3,1);
C{1,1} = A{1};
C{2,1} = b{1};
C{3,1} = c(1);
Plot2DEllipsoid(C,m);
hold on;

plot(x(1),x(2),'r*');
hold on;

plot(x0(1),x0(2),'b*');
hold on;

