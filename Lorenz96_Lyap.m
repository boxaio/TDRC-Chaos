clc,clear;

% dimension 
n = 5;

F = 8;

a = diag(ones(1, n - 1), 1) + diag(-ones(1, n - 2), -2);
a(1, n - 1) = -1;
a(2, n) = -1;
a(n, 1) = 1;
b = diag(ones(1, n - 1), -1);
b(1, n) = 1;

Lorenz96_fun = @(x)(a * x).*(b * x) - x + F;

syms xs [1 n]

g = eval(['@(xs1,xs2,xs3,xs4,xs5)' vectorize(jacobian(Lorenz96_fun(vec(xs)),xs))]);

% Lorenz96_jaco = @(x)g(x(1),x(2),x(3),x(4),x(5));

Lorenz96_grad = @(x, y)vec(g(x(1),x(2),x(3),x(4),x(5)) * y);


dt = 1e-4;
nmax = 4e7; 
norm_steps = 1e4;  % ensure : norm_steps * dt = 1

num_lyaps = n; 

x = [-0.5; 0.2; 1.5; -0.1; 2.3];

y = eye(n);

trans_it = nmax * 0.1;

count = 0;
Rii = zeros(num_lyaps, nmax / norm_steps);

for i = 1 : nmax
    % Evolve Lorenz system and its tangent map forward for 'num_steps' steps
    x = x + dt * Lorenz96_fun(x);
    y = y + dt * reshape(Lorenz96_grad(x, y), n, n);
    % Normalize tangent vectors and record normalization
    % QR decomposition, Q is of size(n, n), R is of size(n, num_lyaps)
    if mod(i, norm_steps) == 0
        count = count + 1;
        % QR decomposition, Q is of size(N, N), R is of size(N, num_lyaps)
        [Q, R] = qr(y, 0);
        Rii(:, count) = log(diag(R));
        % new Y
        y = Q(:, 1:num_lyaps);
    end
end

le = real(sum(Rii, 2)./(nmax * dt)); 

% (0.49213, 0.0, -0.5249, -1.2815, -3.6775)

for j = 1 : n
    disp(['Lambda_',num2str(j),' = ',num2str(le(j))])
end














