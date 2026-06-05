%% ============================================================
%  SOLAR real-data analysis: M2 Levine residual acceleration
%
%  Expected folder structure:
%
%  Project main folder/
%  ├── Current MATLAB script
%  ├── Age methylation data/
%  │     ├── X_matrix_final.mat
%  │     └── Pheno_for_MATLAB_with_epi_ages.csv
%  └── SOLAR supp funs/
%
%  Output:
%      Age methylation data/Output/
%
%  All CSV files begin with:
%      Output_from_M2_....
% ============================================================


% changes specific for real data analysis: Num_iters = 500, W_stop = 250, c_kappa = 1, q_max  = 15


clearvars -except X
clc;
rng(123);

fprintf('MATLAB mem (Windows only):\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Relative directories
% ============================================================

data_dir = 'Age methylation data';
supp_dir = 'SOLAR supp funs';

out_dir = 'Output';

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if ~exist(supp_dir, 'dir')
    error('SOLAR support-function folder not found: %s', supp_dir);
end

addpath(supp_dir);

%% ============================================================
%  File names
% ============================================================

X_mat_file = fullfile(data_dir, 'X_matrix_final.mat');

pheno_file = fullfile( ...
    data_dir, ...
    'Pheno_for_MATLAB_with_epi_ages.csv' ...
);

if ~exist(X_mat_file, 'file')
    error('X matrix file not found: %s', X_mat_file);
end

if ~exist(pheno_file, 'file')
    error('Phenotype file not found: %s', pheno_file);
end

%% ============================================================
%  Load phenotype data
% ============================================================

if exist('pheno_FULL', 'var') && exist('y_raw_FULL', 'var')
    
    fprintf('\nUsing phenotype table and Levine outcome already available in memory.\n');
    
    pheno = pheno_FULL;
    w_raw = y_raw_FULL;
    
else
    
    fprintf('\nReading phenotype file...\n');
    
    pheno_FULL = readtable( ...
        pheno_file, ...
        'VariableNamingRule', 'preserve' ...
    );
    
    if ~ismember( ...
            'AgeAccelResidual_Levine', ...
            pheno_FULL.Properties.VariableNames)
        
        error(['AgeAccelResidual_Levine column ', ...
               'not found in phenotype file.']);
    end
    
    y_raw_FULL = pheno_FULL.AgeAccelResidual_Levine;
    
    pheno = pheno_FULL;
    w_raw = y_raw_FULL;
    
end

%% ============================================================
%  Load methylation matrix X
%  Important: only one full X-sized matrix is kept in memory.
% ============================================================

if exist('X', 'var')
    
    fprintf('\nUsing X already available in memory.\n');
    
else
    
    fprintf('\nLoading X_matrix_final.mat. This may take time...\n');
    
    tic;
    
    tmp = load(X_mat_file);
    var_names = fieldnames(tmp);
    
    if numel(var_names) ~= 1
        
        fprintf('\nVariables found in X_matrix_final.mat:\n');
        disp(var_names);
        
        error(['Expected exactly one numeric matrix ', ...
               'inside X_matrix_final.mat.']);
    end
    
    X = tmp.(var_names{1});
    clear tmp
    
    fprintf('Finished loading X in %.2f seconds.\n', toc);
end

%% ============================================================
%  Basic dimension checks
% ============================================================

[n, p] = size(X);

fprintf('\nLoaded methylation matrix:\n');
fprintf('n = %d samples\n', n);
fprintf('p = %d CpGs\n', p);

fprintf('\nPhenotype table:\n');
fprintf('n_pheno = %d rows\n', height(pheno));

if height(pheno) ~= n
    error(['Number of phenotype rows (%d) does not match ', ...
           'number of X rows (%d).'], ...
           height(pheno), n);
end

if length(w_raw) ~= n
    error(['Length of Levine residual outcome ', ...
           'does not match number of X rows.']);
end

%% ============================================================
%  Remove missing outcome rows, if any
%  Note: this subsets X only if missing outcomes exist.
% ============================================================

valid_idx = ~isnan(w_raw);

if any(~valid_idx)
    
    fprintf(['\nRemoving %d samples with missing ', ...
             'Levine residual acceleration.\n'], ...
             sum(~valid_idx));
    
    X = X(valid_idx, :);
    w_raw = w_raw(valid_idx);
    pheno = pheno(valid_idx, :);
end

[n, p] = size(X);

fprintf('\nAnalysis data after filtering:\n');
fprintf('n = %d samples\n', n);
fprintf('p = %d CpGs\n', p);

%% ============================================================
%  Center X and Levine residual acceleration
%  Memory-conscious column-wise centering of X.
% ============================================================

fprintf('\nCentering X and Levine residual acceleration...\n');

center_tic = tic;

X_mean = mean(X, 1);
w_mean = mean(w_raw);

% Memory-conscious in-place centering.
% This avoids intentionally creating a second X-sized matrix.
for jj = 1:p
    X(:, jj) = X(:, jj) - X_mean(jj);
end

w = w_raw - w_mean;

center_time_sec = toc(center_tic);

fprintf('Finished centering in %.2f seconds.\n', center_time_sec);

Xnorm2 = norm(X, 'fro')^2;
Xnorm  = sqrt(Xnorm2);

%% ============================================================
%  SOLAR settings
% ============================================================

method_name = 'SOLAR_Levine_full';

print_every = 2;

q_init = 10;
q_min  = 1;
q_max  = 15;

Num_iters = 500;

perform_gram_svd = 1;

T0 = 1.0;
T_min = 1e-2;
cool_pow = 0.7;

rho2 = 0.5^2;
g2   = 25.0;

%% ============================================================
%  Rank-penalty tuning
% ============================================================

c_kappa = 1;
zeta = 2/3;

kappa = c_kappa * ...
    (((sqrt(n) + sqrt(p))^2) / (n + p))^zeta;

q_max_eff  = min([q_max, n - 1, p]);
q_init_eff = min(q_init, q_max_eff);

fprintf('\nSOLAR tuning/settings:\n');
fprintf('method_name = %s\n', method_name);
fprintf('q_init = %d | q_min = %d | q_max = %d\n', q_init, q_min, q_max);
fprintf('q_init_eff = %d | q_max_eff = %d\n', q_init_eff, q_max_eff);
fprintf('Num_iters = %d | print_every = %d\n', Num_iters, print_every);
fprintf('perform_gram_svd = %d\n', perform_gram_svd);
fprintf('kappa = %.6f\n', kappa);
fprintf('rho2 = %.6f | g2 = %.6f\n', rho2, g2);

%% ============================================================
%  Precompute SVD(X)
% ============================================================

fprintf('\nPrecomputing leading SVD components...\n');

svd_tic = tic;

if perform_gram_svd == 0
    
    [U_full, S_full, V_full] = svds(X, q_max_eff);
    
else
    
    G = double(X * X');
    
    [U_full, S2] = eigs(G, q_max_eff);
    
    s = sqrt(diag(S2));
    S_full = diag(s);
    
    V_full = X' * U_full;
    
    for k = 1:q_max_eff
        V_full(:, k) = V_full(:, k) / s(k);
    end
end

svd_time_sec = toc(svd_tic);

fprintf('Finished SVD/precomputation in %.2f seconds.\n', svd_time_sec);

%% ============================================================
%  Constants
% ============================================================

fprintf('\nEstimating variance constants...\n');

sigma2 = estimate_sigma2_resid_baseline(X, 2);
tau2   = var(w, 1);

prec_lik   = 1 / sigma2;
prec_prior = 1 / rho2;
prec_post  = prec_lik + prec_prior;

tau_inv2 = 1 / tau2;
g_inv2   = 1 / g2;
prec_b   = tau_inv2 + g_inv2;
beta_shrink = tau_inv2 / prec_b;

fprintf('sigma2_est = %.6f\n', sigma2);
fprintf('tau2_est   = %.6f\n', tau2);

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
%  Best-state tracking
% ============================================================

obj_best  = -Inf;
q_best    = q;
H_best    = H;
V_best    = V;
D_best    = D;
beta_best = beta;

recon_best = Inf;

%% ============================================================
%  Early-stopping parameters
% ============================================================

W_stop      = min(250, floor(0.5 * Num_iters));
epsJ_stop   = 1e-6;
epsAcc_stop = 0.01;
etaT_stop   = 0.10;

best_obj_hist = nan(Num_iters, 1);
best_q_hist   = nan(Num_iters, 1);
q_hist         = nan(Num_iters, 1);
q_accept_hist  = zeros(Num_iters, 1);

%% ============================================================
%  Iterative MAP + trans-dimensional search
% ============================================================

fprintf('\nStarting SOLAR optimization...\n');

opt_tic = tic;
last_iter = Num_iters;

for iter = 1:Num_iters
    
    T = max(T_min, T0 * (iter^(-cool_pow)));
    
    %% ----- V update -----
    
    Fv = X' * (H * D);
    [U1, ~, V1] = svd(Fv, 'econ');
    V = U1 * V1';
    
    %% ----- H update -----
    
    Fh = (w * beta') / tau2 + (X * V * D) / sigma2;
    [U2, ~, V2] = svd(Fh, 'econ');
    H = U2 * V2';
    
    %% ----- beta update -----
    
    beta = beta_shrink * (H' * w);
    
    %% ----- Canonical supervised orientation -----
    
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
    
    %% ----- Re-diagonalize -----
    
    R = H' * X * V;
    [Ur, Sr, Vr] = svd(R, 'econ');
    
    H = H * Ur;
    V = V * Vr;
    
    d = (prec_lik / prec_post) * diag(Sr);
    D = diag(d);
    
    beta = beta_shrink * (H' * w);
    
    %% ----- Objective without forming Xhat -----
    
    XV  = X * V;
    XVD = XV * D;
    
    residX2 = Xnorm2 + norm(D, 'fro')^2 - 2 * trace(H' * XVD);
    
    like_curr = ...
        -0.5 / sigma2 * residX2 ...
        -0.5 / tau2   * norm(w - H * beta)^2 ...
        -0.5 / rho2   * norm(d, 2)^2 ...
        -0.5 / g2     * norm(beta, 2)^2;
    
    dfq = q * (n + p - 2*q) + 2*q;
    pen_curr = kappa * 0.5 * log(n * p) * dfq;
    
    obj_curr = like_curr - pen_curr;
    
    %% ----- Update best state -----
    
    if obj_curr > obj_best
        
        obj_best  = obj_curr;
        q_best    = q;
        H_best    = H;
        V_best    = V;
        D_best    = D;
        beta_best = beta;
        
        recon_best = sqrt(max(residX2, 0)) / Xnorm;
    end
    
    %% ----- Propose rank move -----
    
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
        
        dfq_prop = q_prop * (n + p - 2*q_prop) + 2*q_prop;
        pen_prop = kappa * 0.5 * log(n * p) * dfq_prop;
        
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
    
    %% ----- Print progress -----
    
    if mod(iter, print_every) == 0
        
        XV  = X * V;
        XVD = XV * D;
        
        residX2_iter = Xnorm2 + norm(D, 'fro')^2 - ...
            2 * trace(H' * XVD);
        
        recon_iter = sqrt(max(residX2_iter, 0)) / Xnorm;
        
        fprintf(['Iter %5d | q=%2d | recon=%.4f | obj=%.3e | ', ...
                 'best(q=%2d,recon=%.4f,obj=%.3e)\n'], ...
            iter, q, recon_iter, obj_curr, ...
            q_best, recon_best, obj_best);
    end
    
    best_obj_hist(iter) = obj_best;
    best_q_hist(iter)   = q_best;
    
    %% ----- Early stopping -----
    
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
            
            fprintf(['Early stop at iter=%d: T=%.2e, q_best=%d stable, ', ...
                     'rel_improve=%.2e, acc_rate=%.3f\n'], ...
                iter, T, q_best, rel_improve, acc_rate);
            
            last_iter = iter;
            break;
        end
    end
end

optimization_time_sec = toc(opt_tic);

%% ============================================================
%  Use best MAP state
% ============================================================

q = q_best;
H = H_best;
V = V_best;
D = D_best;
beta = beta_best;
d = diag(D);

fprintf('\nSOLAR optimization complete.\n');
fprintf('Selected rank q_hat = %d\n', q);
fprintf('Optimization time = %.2f seconds\n', optimization_time_sec);

%% ============================================================
%  Fitted values and validation metrics
% ============================================================

w_hat_centered = H * beta;

w_hat_original = w_hat_centered + w_mean;
w_original     = w + w_mean;

rmse_w = sqrt(mean((w_original - w_hat_original).^2));
mae_w  = mean(abs(w_original - w_hat_original));
corr_w = corr(w_hat_original, w_original, 'Rows', 'complete');

ss_res = sum((w_original - w_hat_original).^2);
ss_tot = sum((w_original - mean(w_original)).^2);
R2_w   = 1 - ss_res / ss_tot;

%% ============================================================
%  CpG importance scores
% ============================================================

fprintf('\nComputing CpG importance scores...\n');

factor_weight = (d.^2) .* abs(beta);

importance = (V.^2) * factor_weight;

importance_sum = sum(importance);

if importance_sum > 0
    importance_scaled = importance / importance_sum;
else
    importance_scaled = importance;
end

[~, idx_sorted] = sort(importance, 'descend');

% ============================================================
% Existing top-200 output (keep unchanged)
% ============================================================

topK = min(200, p);
top_idx = idx_sorted(1:topK);

% ============================================================
% Additional enrichment-analysis CpG sets
% (new outputs only; existing outputs unchanged)
% ============================================================

topK_500  = min(500,  p);
topK_2000 = min(2000, p);
topK_5000 = min(5000, p);

top_idx_500  = idx_sorted(1:topK_500);
top_idx_2000 = idx_sorted(1:topK_2000);
top_idx_5000 = idx_sorted(1:topK_5000);

%% ============================================================
%  Factor-strength quantities
% ============================================================

factor_id = (1:q)';

factor_strength_raw = d(:).^2;
factor_strength_prop = factor_strength_raw / sum(factor_strength_raw);

supervised_factor_weight = factor_weight(:);
supervised_factor_weight_prop = ...
    supervised_factor_weight / sum(supervised_factor_weight);

%% ============================================================
%  Prepare minimal output tables
% ============================================================

fprintf('\nPreparing output CSV files...\n');

%% ----- Summary table -----

summary_tbl = table;

summary_tbl.method = string(method_name);
summary_tbl.outcome = string('AgeAccelResidual_Levine');

summary_tbl.n = n;
summary_tbl.p = p;

summary_tbl.q_hat = q;
summary_tbl.q_init = q_init;
summary_tbl.q_min = q_min;
summary_tbl.q_max = q_max;

summary_tbl.Num_iters_requested = Num_iters;
summary_tbl.Num_iters_completed = last_iter;

summary_tbl.rmse_in_sample = rmse_w;
summary_tbl.mae_in_sample  = mae_w;
summary_tbl.corr_in_sample = corr_w;
summary_tbl.R2_in_sample   = R2_w;

summary_tbl.recon_best = recon_best;
summary_tbl.obj_best = obj_best;

summary_tbl.sigma2_est = sigma2;
summary_tbl.tau2_est = tau2;
summary_tbl.kappa = kappa;
summary_tbl.rho2 = rho2;
summary_tbl.g2 = g2;

summary_tbl.center_time_sec = center_time_sec;
summary_tbl.svd_time_sec = svd_time_sec;
summary_tbl.optimization_time_sec = optimization_time_sec;
summary_tbl.total_runtime_sec = ...
    center_time_sec + svd_time_sec + optimization_time_sec;

%% ----- Factor-strength table -----

factor_tbl = table( ...
    factor_id, ...
    d(:), ...
    beta(:), ...
    factor_strength_raw, ...
    factor_strength_prop, ...
    supervised_factor_weight, ...
    supervised_factor_weight_prop, ...
    'VariableNames', { ...
        'factor', ...
        'd', ...
        'beta', ...
        'd_squared', ...
        'd_squared_prop', ...
        'd_squared_abs_beta', ...
        'd_squared_abs_beta_prop' ...
    } ...
);

%% ----- Subject-level latent scores and fitted values -----

score_tbl = table;

score_tbl.sample_index = (1:n)';
score_tbl.Age_months = pheno.Age_months;
score_tbl.y_observed = w_original;
score_tbl.y_fitted = w_hat_original;
score_tbl.residual = w_original - w_hat_original;

for kk = 1:q
    score_tbl.(sprintf('H_factor_%02d', kk)) = H(:, kk);
end

%% ----- Top CpG importance table -----

top_tbl = table;

top_tbl.rank = (1:topK)';
top_tbl.CpG_index = top_idx(:);
top_tbl.importance = importance(top_idx);
top_tbl.importance_scaled = importance_scaled(top_idx);

for kk = 1:q
    top_tbl.(sprintf('loading_V_factor_%02d', kk)) = ...
        V(top_idx, kk);
end

for kk = 1:q
    top_tbl.(sprintf('importance_component_factor_%02d', kk)) = ...
        (V(top_idx, kk).^2) * factor_weight(kk);
end

%% ----- Top 500 CpGs for enrichment analysis -----

top_tbl_500 = table;

top_tbl_500.rank = (1:topK_500)';
top_tbl_500.CpG_index = top_idx_500(:);
top_tbl_500.importance = importance(top_idx_500);
top_tbl_500.importance_scaled = importance_scaled(top_idx_500);

%% ----- Top 2000 CpGs for enrichment analysis -----

top_tbl_2000 = table;

top_tbl_2000.rank = (1:topK_2000)';
top_tbl_2000.CpG_index = top_idx_2000(:);
top_tbl_2000.importance = importance(top_idx_2000);
top_tbl_2000.importance_scaled = importance_scaled(top_idx_2000);

%% ----- Top 5000 CpGs for enrichment analysis -----

top_tbl_5000 = table;

top_tbl_5000.rank = (1:topK_5000)';
top_tbl_5000.CpG_index = top_idx_5000(:);
top_tbl_5000.importance = importance(top_idx_5000);
top_tbl_5000.importance_scaled = importance_scaled(top_idx_5000);

%% ----- Minimal optimization trace -----

trace_tbl = table;

trace_tbl.iter = (1:last_iter)';
trace_tbl.q_current = q_hist(1:last_iter);
trace_tbl.q_best = best_q_hist(1:last_iter);
trace_tbl.obj_best = best_obj_hist(1:last_iter);
trace_tbl.q_accepted = q_accept_hist(1:last_iter);

%% ============================================================
%  Save CSV outputs
% ============================================================

writetable( ...
    summary_tbl, ...
    fullfile(out_dir, 'Output_from_M2_summary.csv') ...
);

writetable( ...
    factor_tbl, ...
    fullfile(out_dir, 'Output_from_M2_factor_strength.csv') ...
);

writetable( ...
    score_tbl, ...
    fullfile(out_dir, 'Output_from_M2_latent_scores_and_fitted.csv') ...
);

writetable( ...
    top_tbl, ...
    fullfile(out_dir, 'Output_from_M2_top_CpG_importance.csv') ...
);

writetable( ...
    top_tbl_500, ...
    fullfile(out_dir, ...
    'Output_from_M2_top_CpG_importance_for_enrichment_500.csv') ...
);

writetable( ...
    top_tbl_2000, ...
    fullfile(out_dir, ...
    'Output_from_M2_top_CpG_importance_for_enrichment_2000.csv') ...
);

writetable( ...
    top_tbl_5000, ...
    fullfile(out_dir, ...
    'Output_from_M2_top_CpG_importance_for_enrichment_5000.csv') ...
);

writetable( ...
    trace_tbl, ...
    fullfile(out_dir, 'Output_from_M2_optimization_trace.csv') ...
);

fprintf('\nSaved CSV outputs in folder: %s\n', out_dir);

fprintf('\nKey results:\n');
fprintf('Selected rank q_hat = %d\n', q);
fprintf('In-sample RMSE = %.4f\n', rmse_w);
fprintf('In-sample correlation = %.4f\n', corr_w);
fprintf('In-sample R2 = %.4f\n', R2_w);
fprintf('Total runtime = %.2f seconds\n', ...
    center_time_sec + svd_time_sec + optimization_time_sec);

try
    disp(memory);
catch
end