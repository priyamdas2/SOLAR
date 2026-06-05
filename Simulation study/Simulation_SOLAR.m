clear all;
rng(123);

fprintf('MATLAB mem (Windows only):\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Add SOLAR support functions
% ============================================================

addpath('SOLAR supp funs');

%% ============================================================
%  Simulation scenario
% ============================================================

n = 500;          % 100, 1000
p = 1000;         % 1000, 10000, 100000

q_true = 3;       % 3, 5

train_frac = 0.80;

target_SNR_X = 0.05;   % 0.01, 0.05, 0.2
target_SNR_w = 0.2;    % 0.05, 0.2, 0.5

rep_vec = 1:10;

% kappa = 0.15;    % Smaller kappa => weaker penalty => larger selected q
c_kappa = 0.1;
zeta = 2/3;
kappa = c_kappa * (((sqrt(train_frac*n) + sqrt(p))^2) / (train_frac*n + p))^zeta;



%% ============================================================
%  SOLAR settings
% ============================================================

method_name = 'SOLAR';

print_every = 200;

q_init = 10;
q_min  = 1;
q_max  = 10;

Num_iters = 2000;

perform_gram_svd = 1;

T0 = 1.0;
T_min = 1e-2;
cool_pow = 0.7;



rho2 = 0.5^2;
g2   = 25.0;


%% ============================================================
%  Directories
% ============================================================

data_dir = 'Data';
out_dir  = 'Output';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ============================================================
%  Filename strings
% ============================================================

snrX_str = strrep(sprintf('%.3f', target_SNR_X), '.', 'p');
snrw_str = strrep(sprintf('%.3f', target_SNR_w), '.', 'p');

%% ============================================================
%  Replication loop
% ============================================================

for rep = rep_vec
    
    fprintf('\n=== %s: Replication %d ===\n', method_name, rep);
    
    %% ============================================================
    %  Load simulated data
    % ============================================================
    
    fname_in = sprintf( ...
        ['Simulated_data_snrX_%s_snrw_%s_qTrue_%d_' ...
        'n_%d_p_%d_rep_%02d.mat'], ...
        snrX_str, snrw_str, q_true, n, p, rep);
    
    data_path = fullfile(data_dir, fname_in);
    
    if ~exist(data_path, 'file')
        warning('File not found: %s. Skipping.', data_path);
        continue;
    end
    
    tmp_load = load(data_path, 'simdata');
    simdata = tmp_load.simdata;
    
    X = simdata.X;
    w = simdata.w;
    
    H_true = simdata.H_true;
    V_true = simdata.V_true;
    D_true = simdata.D_true;
    beta_true = simdata.beta_true;
    
    [n_check, p_check] = size(X);
    
    if n_check ~= n || p_check ~= p
        error('Loaded data dimensions do not match specified n and p.');
    end
    
    %% ============================================================
    %  Train/test split
    % ============================================================
    % Same split rule as PCR/PLS/SPC: rng(1000 + rep), 80/20 split.
    
    rng(1000 + rep);
    
    idx = randperm(n);
    n_train = round(train_frac * n);
    
    train_idx = idx(1:n_train);
    test_idx  = idx(n_train+1:end);
    
    X_train = X(train_idx, :);
    w_train = w(train_idx);
    
    X_test = X(test_idx, :);
    w_test = w(test_idx);
    
    H_true_train = H_true(train_idx, :);
    
    %% ============================================================
    %  Track runtime
    % ============================================================
    
    tic;
    
    %% ============================================================
    %  Center training and test data using training means only
    % ============================================================
    
    X_mean = mean(X_train, 1);
    w_mean = mean(w_train);
    
    X_train_c = X_train - X_mean;
    w_train_c = w_train - w_mean;
    
    X_test_c = X_test - X_mean;
    
    % Use centered training data inside SOLAR
    X = X_train_c;
    w = w_train_c;
    
    [n_fit, p_fit] = size(X);
    
    q_max_eff = min([q_max, n_fit-1, p_fit]);
    q_init_eff = min(q_init, q_max_eff);
    
    Xnorm2 = norm(X, 'fro')^2;
    Xnorm  = sqrt(Xnorm2);
    
    %% ============================================================
    %  Precompute SVD(X)
    % ============================================================
    
    if perform_gram_svd == 0
        
        % --- ORIGINAL METHOD ---
        [U_full, S_full, V_full] = svds(X, q_max_eff);
        
    else
        
        % --- GRAM MATRIX SVD ---
        G = X * X';
        
        [U_full, S2] = eigs(G, q_max_eff);
        
        s = sqrt(diag(S2));
        S_full = diag(s);
        
        V_full = X' * U_full;
        
        for k = 1:q_max_eff
            V_full(:, k) = V_full(:, k) / s(k);
        end
        
    end
    
    %% ============================================================
    %  Constants
    % ============================================================
    
    
    
    sigma2 = estimate_sigma2_resid_baseline(X, 2);
    
    % Crude tau2 initialization from centered training response
    tau2 = var(w, 1);
    
    prec_lik   = 1 / sigma2;
    prec_prior = 1 / rho2;
    prec_post  = prec_lik + prec_prior;
    
    tau_inv2 = 1 / tau2;
    g_inv2   = 1 / g2;
    prec_b   = tau_inv2 + g_inv2;
    beta_shrink = tau_inv2 / prec_b;
    
    %% ============================================================
    %  Initialize
    % ============================================================
    
    q = q_init_eff;
    
    H = U_full(:, 1:q);
    V = V_full(:, 1:q);
    D = (prec_lik / prec_post) * S_full(1:q, 1:q);
    d = diag(D);
    
    beta = beta_shrink * (H' * w);
    
    %% ============================================================
    %  Best-state tracking for MAP
    % ============================================================
    
    obj_best  = -Inf;
    q_best    = q;
    H_best    = H;
    V_best    = V;
    D_best    = D;
    beta_best = beta;
    
    %% ============================================================
    %  Early-stopping parameters
    % ============================================================
    
    W_stop      = 1000;
    epsJ_stop   = 1e-6;
    epsAcc_stop = 0.01;
    etaT_stop   = 0.10;
    
    best_obj_hist = nan(Num_iters, 1);
    best_q_hist   = nan(Num_iters, 1);
    q_hist         = nan(Num_iters, 1);
    q_accept_hist  = zeros(Num_iters, 1);
    
    recon_best = Inf;
    
    %% ============================================================
    %  Iterative MAP + trans-dimensional search
    % ============================================================
    
    for iter = 1:Num_iters
        
        T = max(T_min, T0 * (iter^(-cool_pow)));
        
        %% ----- Within-q alternating MAP -----
        
        % V update
        Fv = X' * (H * D);
        [U1, ~, V1] = svd(Fv, 'econ');
        V = U1 * V1';
        
        % H update
        Fh = (w * beta') / tau2 + (X * V * D) / sigma2;
        [U2, ~, V2] = svd(Fh, 'econ');
        H = U2 * V2';
        
        % beta MAP
        beta = beta_shrink * (H' * w);
        
        % Canonical within-subspace rotation to lock the w direction
        u = H' * w;
        
        if norm(u) > 1e-12
            
            a = u / norm(u);
            e1 = zeros(q, 1);
            e1(1) = 1;
            
            vv = a - e1;
            
            if norm(vv) > 1e-12
                vv = vv / norm(vv);
                Rw = eye(q) - 2 * (vv * vv');
            else
                Rw = eye(q);
            end
            
            H = H * Rw';
            V = V * Rw';
            beta = Rw * beta;
            
        end
        
        % Re-diagonalize
        R = H' * X * V;
        [Ur, Sr, Vr] = svd(R, 'econ');
        
        H = H * Ur;
        V = V * Vr;
        
        d = (prec_lik / prec_post) * diag(Sr);
        D = diag(d);
        
        beta = beta_shrink * (H' * w);
        
        %% ----- Objective without forming full Xhat -----
        
        XV  = X * V;
        XVD = XV * D;
        
        residX2 = Xnorm2 + norm(D, 'fro')^2 - 2 * trace(H' * XVD);
        
        like_curr = ...
            -0.5 / sigma2 * residX2 ...
            -0.5 / tau2   * norm(w - H * beta)^2 ...
            -0.5 / rho2   * norm(d, 2)^2 ...
            -0.5 / g2     * norm(beta, 2)^2;
        
        dfq = q * (n_fit + p_fit - 2*q) + 2*q;
        pen_curr = kappa * 0.5 * log(n_fit * p_fit) * dfq;
        
        obj_curr = like_curr - pen_curr;
        
        %% ----- Update BEST state -----
        
        if obj_curr > obj_best
            
            obj_best  = obj_curr;
            q_best    = q;
            H_best    = H;
            V_best    = V;
            D_best    = D;
            beta_best = beta;
            
            recon_best = sqrt(max(residX2, 0)) / Xnorm;
            
        end
        
        %% ----- Propose q move -----
        
        q_hist(iter) = q;
        
        if rand < 0.5
            
            if q == q_min
                q_prop = q + 1;
            elseif q == q_max_eff
                q_prop = q - 1;
            else
                q_prop = q + (2 * (rand < 0.5) - 1);
            end
            
            H_prop = U_full(:, 1:q_prop);
            V_prop = V_full(:, 1:q_prop);
            D_prop = (prec_lik / prec_post) * S_full(1:q_prop, 1:q_prop);
            d_prop = diag(D_prop);
            
            beta_prop = beta_shrink * (H_prop' * w);
            
            XVp  = X * V_prop;
            XVDp = XVp * D_prop;
            
            residXp2 = Xnorm2 + norm(D_prop, 'fro')^2 - ...
                2 * trace(H_prop' * XVDp);
            
            like_prop = ...
                -0.5 / sigma2 * residXp2 ...
                -0.5 / tau2   * norm(w - H_prop * beta_prop)^2 ...
                -0.5 / rho2   * norm(d_prop, 2)^2 ...
                -0.5 / g2     * norm(beta_prop, 2)^2;
            
            dfq_prop = q_prop * (n_fit + p_fit - 2*q_prop) + 2*q_prop;
            pen_prop = kappa * 0.5 * log(n_fit * p_fit) * dfq_prop;
            
            obj_prop = like_prop - pen_prop;
            
            dObj = obj_prop - obj_curr;
            
            if dObj >= 0 || rand < exp(dObj / T)
                
                q = q_prop;
                H = H_prop;
                V = V_prop;
                D = D_prop;
                d = d_prop;
                beta = beta_prop;
                
                q_accept_hist(iter) = 1;
                q_hist(iter) = q_prop;
                
            end
        end
        
        if mod(iter, print_every) == 0
            
            XV  = X * V;
            XVD = XV * D;
            
            residX2 = Xnorm2 + norm(D, 'fro')^2 - 2 * trace(H' * XVD);
            recon_iter = sqrt(max(residX2, 0)) / Xnorm;
            
            fprintf('Iter %5d | q=%2d | tau2 = %4f | recon=%.4f | obj=%.3e | best(q=%2d,recon=%.4f,obj=%.3e)\n', ...
                iter, q, tau2, recon_iter, obj_curr, ...
                q_best, recon_best, obj_best);
            
        end
        
        best_obj_hist(iter) = obj_best;
        best_q_hist(iter)   = q_best;
        
        %% ============================================================
        %  Early stopping check
        % ============================================================
        
        if iter >= W_stop
            
            at_floor = (T <= (1 + etaT_stop) * T_min);
            
            q_win = best_q_hist(iter-W_stop+1:iter);
            q_stable = all(q_win == q_win(end));
            
            J_now  = best_obj_hist(iter);
            J_prev = best_obj_hist(iter-W_stop+1);
            
            rel_improve = abs(J_now - J_prev) / (1 + abs(J_now));
            J_stable = (rel_improve <= epsJ_stop);
            
            acc_rate = mean(q_accept_hist(iter-W_stop+1:iter));
            acc_small = (acc_rate <= epsAcc_stop);
            
            if at_floor && q_stable && J_stable && acc_small
                
                fprintf('Early stop (SA stage) at iter=%d: T=%.2e, q_best=%d stable, rel_improve=%.2e, acc_rate=%.3f\n', ...
                    iter, T, q_best, rel_improve, acc_rate);
                
                break;
                
            end
        end
    end
    
    %% ============================================================
    %  Use BEST MAP state
    % ============================================================
    
    q = q_best;
    H = H_best;
    V = V_best;
    D = D_best;
    beta = beta_best;
    
    %% ============================================================
    %  Test prediction
    % ============================================================
    
    d_hat = diag(D);
    d_hat_safe = d_hat;
    d_hat_safe(abs(d_hat_safe) < 1e-10) = 1e-10;
    
    H_test = X_test_c * V;
    H_test = H_test ./ (ones(size(H_test, 1), 1) * d_hat_safe');
    
    w_hat_test_centered = H_test * beta;
    w_hat_test = w_hat_test_centered + w_mean;
    
    %% ============================================================
    %  Evaluation metrics
    % ============================================================
    
    q_cmp = min(q, q_true);
    
    H_est_cmp = H(:, 1:q_cmp);
    V_est_cmp = V(:, 1:q_cmp);
    D_est_cmp = D(1:q_cmp, 1:q_cmp);
    
    H_true_cmp = H_true_train(:, 1:q_cmp);
    V_true_cmp = V_true(:, 1:q_cmp);
    D_true_cmp = D_true(1:q_cmp, 1:q_cmp);
    
    %% ----- Rank recovery -----
    
    q_hat = q;
    
    rank_correct = double(q_hat == q_true);
    rank_abs_error = abs(q_hat - q_true);
    
    %% ----- H-subspace recovery on training sample -----
    
    H_projector_error = norm( ...
        H_true_cmp * H_true_cmp' - H_est_cmp * H_est_cmp', ...
        'fro');
    
    %% ----- V-subspace recovery -----
    
    MV = V_true_cmp' * V_est_cmp;
    V_projector_error = sqrt(2*q_cmp - 2*norm(MV, 'fro')^2);
    
    %% ----- Signal reconstruction error for X_signal on training sample -----
    
    inner_X = trace( ...
        D_true_cmp * (H_true_cmp' * H_est_cmp) * ...
        D_est_cmp  * (V_est_cmp' * V_true_cmp));
    
    signal_err2 = norm(D_true_cmp, 'fro')^2 + ...
        norm(D_est_cmp, 'fro')^2 - ...
        2 * inner_X;
    
    X_signal_recon_error = ...
        sqrt(max(signal_err2, 0)) / norm(D_true_cmp, 'fro');
    
    X_signal_rmse = ...
        sqrt(max(signal_err2, 0)) / sqrt(n_train * p);
    
    %% ----- Supervised signal recovery on training sample -----
    
    w_true_signal_train = H_true_train * beta_true;
    w_hat_train_centered = H * beta;
    
    supervised_signal_corr = corr(w_hat_train_centered, w_true_signal_train, ...
        'Rows', 'complete');
    
    supervised_signal_rel_error = ...
        norm(w_hat_train_centered - w_true_signal_train) / norm(w_true_signal_train);
    
    supervised_signal_rmse = ...
        norm(w_hat_train_centered - w_true_signal_train) / sqrt(n_train);
    
    %% ----- Out-of-sample prediction performance on test sample -----
    
    rmse_w = sqrt(mean((w_test - w_hat_test).^2));
    mae_w = mean(abs(w_test - w_hat_test));
    corr_w = corr(w_hat_test, w_test, 'Rows', 'complete');
    
    ss_res = sum((w_test - w_hat_test).^2);
    ss_tot = sum((w_test - mean(w_test)).^2);
    
    R2_w = 1 - ss_res / ss_tot;
    
    %% ----- Beta recovery after Procrustes alignment -----
    
    beta_l2_error = NaN;
    
    if q == q_true
        
        [Uo, ~, Vo] = svd(H_true_train' * H, 'econ');
        R_align = Uo * Vo';
        
        beta_aligned = R_align' * beta;
        
        beta_l2_error = norm(beta_aligned - beta_true);
        
    end
    
    %% ----- Timing -----
    
    runtime_sec = toc;
    
    %% ============================================================
    %  Save one-row CSV output
    % ============================================================
    
    result = table;
    
    result.method = string(method_name);
    
    result.q_true = q_true;
    result.q_hat  = q_hat;
    
    result.rank_correct = rank_correct;
    result.rank_abs_error = rank_abs_error;
    
    result.H_projector_error = H_projector_error;
    result.V_projector_error = V_projector_error;
    
    result.X_signal_recon_error = X_signal_recon_error;
    result.X_signal_rmse = X_signal_rmse;
    
    result.supervised_signal_rmse = supervised_signal_rmse;
    result.supervised_signal_corr = supervised_signal_corr;
    result.supervised_signal_rel_error = supervised_signal_rel_error;
    
    result.rmse_w = rmse_w;
    result.mae_w = mae_w;
    result.corr_w = corr_w;
    result.R2_w = R2_w;
    
    result.beta_l2_error = beta_l2_error;
    
    result.runtime_sec = runtime_sec;
    
    out_fname = sprintf( ...
        ['Output_%s_snrX_%s_snrw_%s_qTrue_%d_' ...
        'n_%d_p_%d_rep_%02d.csv'], ...
        method_name, snrX_str, snrw_str, q_true, n, p, rep);
    
    writetable(result, fullfile(out_dir, out_fname));
    
    fprintf('Saved: %s\n', out_fname);
    fprintf('q_hat = %d | Test RMSE(w) = %.4f | Signal corr = %.4f | Runtime = %.2f sec\n', ...
        q_hat, rmse_w, supervised_signal_corr, runtime_sec);
    
end