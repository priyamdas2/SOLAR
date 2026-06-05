clear all;
rng(123);

fprintf('MATLAB mem (Windows only):\n');
try
    disp(memory);
catch
end

%% ============================================================
%  Simulation scenario
% ============================================================

n = 500;          % 100, 1000
p = 1000;         % 1000, 10000, 100000

q_true = 5;       % 3, 5

target_SNR_X = 0.05;   % 0.01, 0.05, 0.2
target_SNR_w = 0.2;    % 0.05, 0.2, 0.5

rep_vec = 1:10;

%% ============================================================
%  PLS settings
% ============================================================

method_name = 'PLS';

q_grid = 1:10;

train_frac = 0.80;
Kfold = 5;

% Ensure q_grid does not exceed training sample-size limits
n_train_nominal = round(train_frac * n);
q_grid = q_grid(q_grid < n_train_nominal);

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

    loaded_data = load(data_path, 'simdata');
    simdata = loaded_data.simdata;

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

    %% ============================================================
    %  Select number of PLS components by K-fold CV on training data
    % ============================================================

    q_grid_eff = q_grid(q_grid < n_train);

    cv = cvpartition(n_train, 'KFold', Kfold);

    cv_rmse = nan(length(q_grid_eff), 1);

    for ii = 1:length(q_grid_eff)

        q = q_grid_eff(ii);
        fold_mse = nan(Kfold, 1);

        for kk = 1:Kfold

            idx_tr = training(cv, kk);
            idx_va = test(cv, kk);

            X_cv_tr = X_train(idx_tr, :);
            w_cv_tr = w_train(idx_tr);

            X_cv_va = X_train(idx_va, :);
            w_cv_va = w_train(idx_va);

            % Center each CV training fold separately
            X_cv_mean = mean(X_cv_tr, 1);
            w_cv_mean = mean(w_cv_tr);

            X_cv_tr_c = X_cv_tr - X_cv_mean;
            w_cv_tr_c = w_cv_tr - w_cv_mean;

            X_cv_va_c = X_cv_va - X_cv_mean;

            % Fit PLS on CV training fold
            [~, ~, ~, ~, BETA_cv] = plsregress(X_cv_tr_c, w_cv_tr_c, q);

            % Predict validation fold
            w_cv_pred_c = [ones(sum(idx_va), 1), X_cv_va_c] * BETA_cv;
            w_cv_pred = w_cv_pred_c + w_cv_mean;

            fold_mse(kk) = mean((w_cv_va - w_cv_pred).^2);

        end

        cv_rmse(ii) = sqrt(mean(fold_mse, 'omitnan'));

    end

    [~, best_idx] = min(cv_rmse);
    q_hat = q_grid_eff(best_idx);

    %% ============================================================
    %  Final PLS fit on full training data using selected q_hat
    % ============================================================

    [XL_hat, YL_hat, XS_hat, YS_hat, BETA_hat, PCTVAR_hat, MSE_hat, stats_hat] = ...
        plsregress(X_train_c, w_train_c, q_hat);

    % Training fitted values, centered and original scale
    w_hat_train_centered = [ones(n_train, 1), X_train_c] * BETA_hat;
    w_hat_train = w_hat_train_centered + w_mean;

    % Test predictions
    w_hat_test_centered = [ones(length(test_idx), 1), X_test_c] * BETA_hat;
    w_hat_test = w_hat_test_centered + w_mean;

    %% ============================================================
    %  Construct orthonormal PLS subspaces from training fit
    % ============================================================
    % PLS scores/loadings are not the same as an orthogonal SVD factorization.
    % For fair rotation-invariant comparison, we evaluate the column spaces
    % spanned by the PLS scores and X-loadings.

    H_hat = orth(XS_hat(:, 1:q_hat));
    V_hat = orth(XL_hat(:, 1:q_hat));

    % If numerical rank drops, adjust effective rank for stable evaluation
    q_eff = min(size(H_hat, 2), size(V_hat, 2));

    H_hat = H_hat(:, 1:q_eff);
    V_hat = V_hat(:, 1:q_eff);

    % Core matrix for projected X reconstruction on PLS subspaces
    C_hat = H_hat' * X_train_c * V_hat;

    % H-scale supervised coefficient
    beta_hat = H_hat' * w_hat_train_centered;

    %% ============================================================
    %  Evaluation metrics
    % ============================================================

    q_cmp = min(q_eff, q_true);

    H_est_cmp = H_hat(:, 1:q_cmp);
    V_est_cmp = V_hat(:, 1:q_cmp);
    C_est_cmp = C_hat(1:q_cmp, 1:q_cmp);

    H_true_cmp = H_true_train(:, 1:q_cmp);
    V_true_cmp = V_true(:, 1:q_cmp);
    D_true_cmp = D_true(1:q_cmp, 1:q_cmp);

    %% ----- Rank recovery -----

    rank_correct  = double(q_hat == q_true);
    rank_abs_error = abs(q_hat - q_true);

    %% ----- H-subspace recovery -----

    H_projector_error = norm( ...
        H_true_cmp * H_true_cmp' - H_est_cmp * H_est_cmp', ...
        'fro');

    %% ----- V-subspace recovery -----

    MV = V_true_cmp' * V_est_cmp;
    V_projector_error = sqrt(2*q_cmp - 2*norm(MV, 'fro')^2);

    %% ----- Signal reconstruction error for X_signal on training sample -----
    % PLS does not directly estimate a diagonal SVD-like D.
    % We therefore reconstruct the best projected low-rank signal over the
    % estimated PLS left/right subspaces using C_hat = H_hat' X_train_c V_hat.
    %
    % Computes:
    % ||H_true_train D_true V_true' - H_hat C_hat V_hat'||_F
    % -----------------------------------------------------------
    %        ||H_true_train D_true V_true'||_F

    inner_X = trace( ...
        D_true_cmp * (H_true_cmp' * H_est_cmp) * ...
        C_est_cmp  * (V_est_cmp' * V_true_cmp));

    signal_err2 = norm(D_true_cmp, 'fro')^2 + ...
        norm(C_est_cmp, 'fro')^2 - ...
        2 * inner_X;
    
    X_signal_recon_error = ...
        sqrt(max(signal_err2, 0)) / norm(D_true_cmp, 'fro');
    
    % Entrywise/root-mean-square reconstruction error
    X_signal_rmse = ...
        sqrt(max(signal_err2, 0)) / sqrt(n_train * p);
    
    %% ----- Supervised signal recovery on training sample -----
    
    w_true_signal_train = H_true_train * beta_true;
    
    supervised_signal_corr = corr(w_hat_train_centered, w_true_signal_train, ...
        'Rows', 'complete');
    
    supervised_signal_rel_error = ...
        norm(w_hat_train_centered - w_true_signal_train) / norm(w_true_signal_train);
    
    % Root-mean-square supervised signal error
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
    % Beta is only directly comparable when q_hat equals q_true and the
    % effective numerical rank also matches q_true.

    beta_l2_error = NaN;

    if q_hat == q_true && q_eff == q_true

        [Uo, ~, Vo] = svd(H_true_train' * H_hat, 'econ');
        R_align = Uo * Vo';

        beta_aligned = R_align' * beta_hat;

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