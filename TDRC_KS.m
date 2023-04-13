function [td, tt, xx, u_true, prediction, varargout] ...
            = TDRC_KS(m, set_average_degree, approx_reservoir_size,...
                      rho, num_delay, gama, beta, c, num_inter_step, ...
                      compute_lyap, num_lyap)


mu = m.mu;
data = m.ut;
L = m.L;
dt = m.dt;

[~, input_dim] = size(data);
xd = L * (-input_dim/2 + 1 : input_dim/2)' / input_dim;

%% TDRC settings

% structure settings
sparsity = set_average_degree / approx_reservoir_size;
nodes_per_input = round(approx_reservoir_size / input_dim);
N = input_dim * nodes_per_input;

W = adj_mat_uniform(N, rho, sparsity);

% store the reservoir state
M = N * num_delay;
xr = zeros(M, 1);

transient_len = 10000;   % discard initial transient
train_len = 40000;       % training length
pred_len = 3000;        % prediction length (for autonomous prediction)

w_in = zeros(N, input_dim); 

% input coupling mask
q = nodes_per_input;
for i = 1 : input_dim
    rng(i)
    ip = gama * (-1 + 2 * rand(q,1));
    w_in((i - 1) * q + 1 : i * q, i) = ip;
end

% activation function
sigma = @(x)tanh(x);

disp(['washout initial transients...'])

for i = 1 : transient_len
    rt = sigma(W * xr(1 : num_delay : N * num_delay) + w_in * data(i,:)');
    xr_temp = reshape(xr, N, num_delay);
    xr_temp = [rt, xr_temp(:, 1 : end-1)];
    xr = reshape(xr_temp, M, 1);
end

disp(['training...'])


xr_train = zeros(M, train_len);
xr_train(:, 1) = xr;

for i = 1 : train_len - 1
    rt = tanh(W * xr_train(1 : num_delay :M, i) + w_in * data(transient_len+i,:)');
    xr_temp = reshape(xr_train(:, i), N, num_delay);
    xr_temp = [rt, xr_temp(:, 1 : end-1)];
    xr_train(:, i+1) = reshape(xr_temp, M, 1);
end

% train, wout is of size (num_inputs, M)
Wout = data(transient_len + 1 : transient_len + train_len, :)' ...
          * xr_train' * pinv(xr_train * xr_train' + beta * eye(M));

err_train = norm(data(transient_len+1:transient_len+train_len,:)' - Wout * xr_train) / train_len;
disp(['Training error : ', num2str(err_train)])

disp(['autonomous prediction with intermittent coupling to the source system...'])

xr_pred = xr_train(:, end);

prediction = zeros(pred_len, input_dim);
prediction(1,:) = Wout * xr_pred;

% initial condition
u0 = prediction(1,:)';
u_true = KS_Solver(L, input_dim, dt, pred_len * dt, mu, u0);

for i = 1 : pred_len - 1
    % obtain true state
    if mod(i, num_inter_step) == 0
        input_to_res = prediction(i, :)' + c * (u_true(i, :)' - prediction(i, :)');
    else
        input_to_res = prediction(i, :)';
    end
    % update reservoir state
    rt = sigma(W * xr_pred(1 : num_delay : M) + w_in * input_to_res);
    xr_temp = reshape(xr_pred, N, num_delay);
    xr_temp = [rt, xr_temp(:, 1 : end-1)];
    xr_pred = reshape(xr_temp, M, 1);
    % reservoir output
    prediction(i+1, :) = Wout * xr_pred;
end

ks_spect = m.spectrum;

td = dt * (1 : pred_len) * ks_spect(1);
[xx, tt] = meshgrid(xd+L/2, td);


if compute_lyap

    %% now compute the Lyapunov exponents of TDRC with intermittent 
    %% access to actual system data

    % maximum number of steps
    nmax = 2e3;
    trans_it = 0.3*nmax;
    % number of steps, greater than num_inter_step
    num_step = 2 * num_inter_step;

    xr_ly = xr_pred;
    % initial condition
    u_0 = Wout * xr_ly;
    u_ly = KS_Solver(L, input_dim, dt, nmax * num_step * dt, mu, u_0);

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
