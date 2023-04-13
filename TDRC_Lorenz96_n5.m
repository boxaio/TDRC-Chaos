function [t_pred, u_target, prediction, varargout] ...
           = TDRC_Lorenz96_n5(dt, set_average_degree, approx_reservoir_size,...
                              rho, num_delay, gama, beta, c, num_inter_step, ...
                              compute_lyap, num_lyap)

n = 5;

F = 8;

a = diag(ones(1, n - 1), 1) + diag(-ones(1, n - 2), -2);
a(1, n - 1) = -1;
a(2, n) = -1;
a(n, 1) = 1;
b = diag(ones(1, n - 1), -1);
b(1, n) = 1;

Lorenz96_fun = @(x)(a * x).*(b * x) - x + F;

Ndata = 50000;

tspan = 0 : dt : (Ndata - 1) * dt;
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[~, x] = ode45(@(t,x)Lorenz96_fun(x), tspan, [-0.5; 0.2; 1.5; -0.1; 2.3], opts);

u = x';

%% TDRC settings

% structure settings
sparsity = set_average_degree / approx_reservoir_size;
nodes_per_input = round(approx_reservoir_size / n);
N = n * nodes_per_input;

W = adj_mat_uniform(N, rho, sparsity);


% store the reservoir state
M = N * num_delay;
xr = zeros(M, 1);

transient_len = 10000;   % discard initial transient
train_len = 40000;       % training length
pred_len = 1000;        % prediction length (for autonomous prediction)

w_in = zeros(N, n); 

% input coupling mask
q = nodes_per_input;
for i = 1 : n
    rng(i)
    ip = gama * (-1 + 2 * rand(q,1));
    w_in((i - 1) * q + 1 : i * q, i) = ip;
end

% activation function
sigma = @(x)tanh(x);

disp(['washout initial transients...'])

for i = 1 : transient_len
    rt = sigma(W * xr(1 : num_delay : N * num_delay) + w_in * u(:, i));
    xr_temp = reshape(xr, N, num_delay);
    xr_temp = [rt, xr_temp(:, 1 : end-1)];
    xr = reshape(xr_temp, M, 1);
end

disp(['training...'])


xr_train = zeros(M, train_len);
xr_train(:, 1) = xr;

for i = 1 : train_len - 1
    rt = sigma(W * xr_train(1 : num_delay : N * num_delay, i)...
               + w_in * u(:, transient_len+i));
    xr_temp = reshape(xr_train(:, i), N, num_delay);
    xr_temp = [rt, xr_temp(:, 1 : end-1)];
    xr_train(:, i+1) = reshape(xr_temp, M, 1);
end

% train, wout is of size (num_inputs, M)
Wout = u(:, transient_len + 1 : transient_len + train_len) ...
          * xr_train' * pinv(xr_train * xr_train' + beta * eye(M));

err_train = norm(u(:,transient_len+1:transient_len+train_len) - Wout * xr_train) / train_len;
disp(['Training error : ', num2str(err_train)])

disp(['autonomous prediction...'])

xr_pred = xr_train(:, end);

prediction = zeros(n, pred_len);
prediction(:, 1) = Wout * xr_pred;

tt = 0 : dt : (pred_len - 1) * dt;
% initial condition
u0 = prediction(:, 1);
[t, ut] = ode45(@(t,x)Lorenz96_fun(x), tt, u0, opts);
u_target = ut';

for i = 1 : pred_len - 1
    % obtain true state
    if mod(i, num_inter_step) == 0
        input_to_res = prediction(:, i) + c * (u_target(:,i) - prediction(:, i));
    else
        input_to_res = prediction(:, i);
    end
    % update reservoir state
    rt = sigma(W * xr_pred(1 : num_delay : M) + w_in * input_to_res);
    xr_temp = reshape(xr_pred, N, num_delay);
    xr_temp = [rt, xr_temp(:, 1 : end-1)];
    xr_pred = reshape(xr_temp, M, 1);
    % reservoir output
    prediction(:, i+1) = Wout * xr_pred;
end

% the largest Lyapunov exponent of Lorenz96
le1 = 0.49213;
t_pred = t * le1;

if compute_lyap

    %% now compute the Lyapunov exponents of TDRC with intermittent 
    %% access to actual system data

    % maximum number of steps
    nmax = 1e3;
    trans_it = 0.3*nmax;
    % number of steps, greater than num_inter_step
    num_step = 2 * num_inter_step;

    xr_ly = xr_pred;
    % initial condition
    u_0 = Wout * xr_ly;
    [~, u_ly] = ode45(@(t,x)Lorenz96_fun(x), 0 : dt : nmax * num_step * dt, u_0, opts);

    out = u_0;

    % initial perturbation vectors (orthogonal), each vector is of length M
    w = orth(rand(M, num_lyap));

    % companion matrix, size (M, M)
    C = zeros(M, M);
    for i = 1 : num_delay - 1
        C(i*N+1:(i+1)*N, (i-1)*N+1:i*N) = eye(N);
    end

    B = [zeros(N, M-N), eye(N)];

    h = zeros(1, num_lyap);
    ly = -1;
    ly_history = [];

    for i = 1 : nmax
        % Evolve reservoir and its tangent map forward for 'num_step' steps
        for j = 1 : num_step
            % obtain true state
            if mod(j, num_inter_step) == 0
                input_to_res = out + c * (u_ly((i-1)*num_step+j, :)' - out);
                W2 = (1 - c) * Wout;
            else
                input_to_res = out;
                W2 = Wout;
            end
            % update reservoir state
            rt = sigma(W * xr_ly(1 : num_delay : M) + w_in * input_to_res);
            xr_temp = reshape(xr_ly, N, num_delay);
            xr_temp = [rt, xr_temp(:, 1 : end-1)];
            xr_ly = reshape(xr_temp, M, 1);
            % reservoir output
            out = Wout * xr_ly;
            % update perturbation vector
            w = C * w + [(1-rt.^2).*((W * B + w_in * W2) * w); zeros(M-N,num_lyap)];
        end
        w_next = w;
        % Normalize tangent vectors and record normalization
        % QR decomposition, Q is of size(n, n), R is of size(n, num_lyaps)
        [Q, R] = qr(w, 0);
        w_next = Q * diag(diag(R));
        h_next = h + log(vecnorm(w_next));
        ly_next = h_next / (i * num_step * dt);
        ly_history = [ly_history; ly_next];
        i
        if i > trans_it && norm(ly_next - ly) < 1e-4
            Le = sort(ly_next, 'descend');
            break;
        end
        h = h_next;
        w = w_next./vecnorm(w_next);
        ly = ly_next;
    end
    
    varargout{1} = ly_next;
    varargout{2} = ly_history;
    
end





















