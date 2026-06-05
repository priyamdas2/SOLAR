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

n = 1000;          % e.g., 1000
p = 10^7;       % e.g., 1000000, 10000000

q_true = 5;        % 3, 5

target_SNR_X = 0.01;   % 0.01, 0.05, 0.2
target_SNR_w = 0.5;    % 0.05, 0.2, 0.5

rep_vec = 1:10;

%% ============================================================
%  Data-generation parameters
% ============================================================

sigma2_true = 0.05^2;

[d_true, info] = generate_diag_from_snr_detectable( ...
    target_SNR_X, n, p, sigma2_true, q_true, 1.15, 1.5);

D_true = diag(d_true);

beta_true = 5 * randn(q_true, 1);

tau2_true = (norm(beta_true, 2)^2) / (n * target_SNR_w);

%% ============================================================
%  SOLAR settings
% ============================================================

method_name = 'SOLAR';

print_every = 2;

q_init = 10;
q_min  = 1;
q_max  = 10;

Num_iters = 400;

perform_gram_svd = 1;

T0 = 1.0;
T_min = 1e-2;
cool_pow = 0.7;

rho2 = 0.5^2;
g2   = 25.0;

%% ============================================================
%  Dimension-adaptive kappa
% ============================================================

c_kappa = 0.1;
zeta = 2/3;

kappa = c_kappa * (((sqrt(n) + sqrt(p))^2) / (n + p))^zeta;

fprintf('Using kappa(n,p) = %.6f\n', kappa);

%% ============================================================
%  Directories
% ============================================================

out_dir = 'Output';

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
    %  Generate simulated data inside replication
    % ============================================================
    
    rng(123 + rep);
    
    [H_true, ~] = qr(randn(n, q_true), 0);
    [V_true, ~] = qr(randn(p, q_true), 0);
    
    % Avoid storing X_signal separately for large p.
    X = H_true * D_true * V_true' + sqrt(sigma2_true) * randn(n, p);
    
    Xnorm2_raw = norm(X, 'fro')^2;
    Xnorm_raw  = sqrt(Xnorm2_raw);
    
    w_signal = H_true * beta_true;
    w = w_signal + sqrt(tau2_true) * randn(n, 1);
    
    % Data-generation diagnostics without forming X_signal explicitly
    snrX_rep = norm(D_true, 'fro')^2 / (n * p * sigma2_true);
    snrw_rep = norm(w_signal, 2)^2 / (n * tau2_true);
    
    fprintf('Target SNR_X = %.4f | Designed SNR_X = %.4f\n', ...
        target_SNR_X, snrX_rep);
    
    fprintf('Target SNR_w = %.4f | Designed SNR_w = %.4f\n', ...
        target_SNR_w, snrw_rep);
    
    %% ============================================================
    %  Track runtime
    % ============================================================
    
    tic;
    
    %% ============================================================
    %  Center full data
    % ============================================================
    
    X_mean = mean(X, 1);
    w_mean = mean(w);
    
    X = X - X_mean;
    w = w - w_mean;
    
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
    
    W_stop      = 200;
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
            
            fprintf('Iter %5d | q=%2d | tau2 = %.4f | recon=%.4f | obj=%.3e | best(q=%2d,recon=%.4f,obj=%.3e)\n', ...
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
    %  Evaluation metrics
    % ============================================================
    
    q_cmp = min(q, q_true);
    
    H_est_cmp = H(:, 1:q_cmp);
    V_est_cmp = V(:, 1:q_cmp);
    D_est_cmp = D(1:q_cmp, 1:q_cmp);
    
    H_true_cmp = H_true(:, 1:q_cmp);
    V_true_cmp = V_true(:, 1:q_cmp);
    D_true_cmp = D_true(1:q_cmp, 1:q_cmp);
    
    %% ----- Rank recovery -----
    
    q_hat = q;
    
    rank_correct = double(q_hat == q_true);
    rank_abs_error = abs(q_hat - q_true);
    
    %% ----- H-subspace recovery -----
    
    H_projector_error = norm( ...
        H_true_cmp * H_true_cmp' - H_est_cmp * H_est_cmp', ...
        'fro');
    
    %% ----- V-subspace recovery -----
    
    MV = V_true_cmp' * V_est_cmp;
    V_projector_error = sqrt(2*q_cmp - 2*norm(MV, 'fro')^2);
    
    %% ----- Signal reconstruction error for X_signal -----
    
    inner_X = trace( ...
        D_true_cmp * (H_true_cmp' * H_est_cmp) * ...
        D_est_cmp  * (V_est_cmp' * V_true_cmp));
    
    signal_err2 = norm(D_true_cmp, 'fro')^2 + ...
        norm(D_est_cmp, 'fro')^2 - ...
        2 * inner_X;
    
    X_signal_recon_error = ...
        sqrt(max(signal_err2, 0)) / norm(D_true_cmp, 'fro');
    
    X_signal_rmse = ...
        sqrt(max(signal_err2, 0)) / sqrt(n * p);
    
    %% ----- Supervised signal recovery -----
    
    w_true_signal = H_true * beta_true;
    w_hat_centered = H * beta;
    
    supervised_signal_corr = corr(w_hat_centered, w_true_signal, ...
        'Rows', 'complete');
    
    supervised_signal_rel_error = ...
        norm(w_hat_centered - w_true_signal) / norm(w_true_signal);
    
    supervised_signal_rmse = ...
        norm(w_hat_centered - w_true_signal) / sqrt(n);
    
    %% ----- In-sample response reconstruction metrics -----
    % These are not out-of-sample prediction metrics.
    
    w_hat = w_hat_centered + w_mean;
    w_original = w + w_mean;
    
    rmse_w = sqrt(mean((w_original - w_hat).^2));
    mae_w = mean(abs(w_original - w_hat));
    corr_w = corr(w_hat, w_original, 'Rows', 'complete');
    
    ss_res = sum((w_original - w_hat).^2);
    ss_tot = sum((w_original - mean(w_original)).^2);
    R2_w = 1 - ss_res / ss_tot;
    
    %% ----- Beta recovery after Procrustes alignment -----
    
    beta_l2_error = NaN;
    
    if q == q_true
        
        [Uo, ~, Vo] = svd(H_true' * H, 'econ');
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
    
    result.snrX_rep = snrX_rep;
    result.snrw_rep = snrw_rep;
    result.kappa = kappa;
    result.sigma2_est = sigma2;
    result.tau2_est = tau2;
    
    out_fname = sprintf( ...
        ['Output_%s_snrX_%s_snrw_%s_qTrue_%d_' ...
        'n_%d_p_%d_rep_%02d.csv'], ...
        method_name, snrX_str, snrw_str, q_true, n, p, rep);
    
    writetable(result, fullfile(out_dir, out_fname));
    
    fprintf('Saved: %s\n', out_fname);
    fprintf('q_hat = %d | In-sample RMSE(w) = %.4f | Signal corr = %.4f | Runtime = %.2f sec\n', ...
        q_hat, rmse_w, supervised_signal_corr, runtime_sec);
    
    %% ============================================================
    %  Clear large objects before next replication
    % ============================================================
    
    clear X w H_true V_true w_signal
    clear U_full S_full V_full G S2 s
    clear H V D beta H_best V_best D_best beta_best
    clear H_est_cmp V_est_cmp D_est_cmp H_true_cmp V_true_cmp D_true_cmp
    clear XV XVD XVp XVDp Fv Fh R Ur Sr Vr U1 V1 U2 V2
    clear X_mean w_mean w_hat w_original w_hat_centered w_true_signal
    clear best_obj_hist best_q_hist q_hist q_accept_hist
    
    try
        disp(memory);
    catch
    end
    
end

%% ============================================================
%  Helper function
% ============================================================

function [d_true, info] = generate_diag_from_snr_detectable(target_SNR_X, n, p, sigma2, q_true, alpha_min, alpha_max)
% Generate singular values d_true with:
%  (i) detectable "shape" relative to noise bulk edge
% (ii) exact Frobenius SNR_X matching target_SNR_X

sigma = sqrt(sigma2);

% Noise bulk edge
bulk_edge = sigma * (sqrt(n) + sqrt(p));

if nargin < 6 || isempty(alpha_min)
    alpha_min = 1.15;
end

if nargin < 7 || isempty(alpha_max)
    alpha_max = alpha_min + 0.35;
end

alpha = linspace(alpha_max, alpha_min, q_true)';
s = bulk_edge * alpha;

S2 = sum(s.^2);
c  = sqrt(target_SNR_X * (n * p * sigma2) / S2);

d_true = c * s;

info.bulk_edge = bulk_edge;
info.alpha     = alpha;
info.scale_c   = c;
info.SNR_check = sum(d_true.^2) / (n * p * sigma2);
info.d_min     = min(d_true);
info.d_max     = max(d_true);

end