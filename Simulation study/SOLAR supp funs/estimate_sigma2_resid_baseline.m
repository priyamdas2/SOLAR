function sigma2_hat = estimate_sigma2_resid_baseline(X, q0)
% Estimate sigma^2 from a baseline rank-q0 approximation of X.
% Uses SVD-based rank-q0 fit and a df correction.

    [n,p] = size(X);
    r = min(n,p);
    if q0 < 0 || q0 >= r
        error('q0 must satisfy 0 <= q0 < min(n,p).');
    end
    if q0 == 0
        % no signal fit: sigma2 from total energy
        sigma2_hat = mean(X(:).^2);
        return;
    end

    [U,S,V] = svd(X,'econ');
    Xhat0 = U(:,1:q0) * S(1:q0,1:q0) * V(:,1:q0)';

    % df for low-rank X part (H,V on Stiefel + diag(D)):
    dfX = q0*(n + p - 2*q0) + q0;

    denom = max(1, n*p - dfX);
    sigma2_hat = norm(X - Xhat0,'fro')^2 / denom;
end