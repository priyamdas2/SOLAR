function tau2_hat = estimate_tau2_resid(w, H, beta)
% tau^2 from regression residuals with df correction.
    n = length(w);
    q = size(H,2);
    denom = max(1, n - q);
    tau2_hat = norm(w - H*beta,2)^2 / denom;
end