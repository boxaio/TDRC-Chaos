clc,clear;

n = 3;

Lorenz63 = @(sigma, rho, beta, x) ...
              [ sigma * (x(2) - x(1));
                x(1) * (rho - x(3)) - x(2);
                x(1) * x(2) - beta * x(3)];

Lorenz63_grad = @(sigma, rho, beta, x, y) ...
                  vec([-sigma, sigma, 0;
                       rho - x(3), -1, -x(1);
                       x(2), x(1), -beta] * y);               
               
sigma = 10;
rho = 28;
beta = 8/3;

dt = 1e-4;
nmax = 4e7; 
norm_steps = 1e4;

num_lyaps = n; 

x = [-0.5,  0.2,  1.5];

y = eye(n);

trans_it = nmax * 0.1;

count = 0;
Rii = zeros(num_lyaps, nmax / norm_steps);

for i = 1 : nmax
    % Evolve Lorenz system and its tangent map forward for 'num_steps' steps
    x = x + dt * Lorenz63(sigma, rho, beta, x);
    y = y + dt * reshape(Lorenz63_grad(sigma, rho, beta, x, y), n, n);
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

% ()


for j = 1 : n
    disp(['Lambda_',num2str(j),' = ',num2str(le(j))])
end









