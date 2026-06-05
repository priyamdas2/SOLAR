function out = compute_kappa_auto_nq(X, sigma2, n, p, tail_frac, q0_for_delta)
% compute_kappa_auto_nq
% q-free kappa_auto using:
%   s_noise2 = median of tail singular values squared
%   Delta_df = df(q0+1)-df(q0) with a fixed baseline q0 (default q0=0)
%
% Inputs:
%   X            : n x p data matrix
%   sigma2       : noise variance sigma^2 (known or estimated)
%   n, p         : dimensions
%   tail_frac    : fraction in (0,1); start index for tail singular values (default 0.5)
%   q0_for_delta : baseline q0 for Delta_df (default 0)
%
% Output struct fields:
%   out.kappa_auto
%   out.s_noise2
%   out.Delta_df
%   out.q0_for_delta
%   out.tail_frac
%   out.tail_idx
%   out.singular_values

    if nargin < 5 || isempty(tail_frac)
        tail_frac = 0.5;
    end
    if nargin < 6 || isempty(q0_for_delta)
        q0_for_delta = 0; % baseline for first-factor increment
    end
    if tail_frac <= 0 || tail_frac >= 1
        error('tail_frac must be in (0,1).');
    end
    if q0_for_delta < 0
        error('q0_for_delta must be >= 0.');
    end
    if q0_for_delta >= min(n,p)
        error('q0_for_delta must be < min(n,p).');
    end

    % --- singular values (economy SVD)
    s = svd(X,'econ');
    r = length(s);

    % --- tail for bulk noise scale
    i0 = max(1, ceil(tail_frac * r));
    tail_idx = i0:r;

    s_noise2 = median(s(tail_idx).^2);

    % --- Delta_df at baseline q0:
    % df(q) = q(n+p-2q)+2q, so Delta_df(q)= (n+p+2)-4q
    Delta_df = (n + p + 2) - 4*q0_for_delta;

    kappa_auto = (s_noise2 / sigma2) / (log(n*p) * Delta_df);

    out = struct();
    out.kappa_auto      = kappa_auto;
    out.s_noise2        = s_noise2;
    out.Delta_df        = Delta_df;
    out.q0_for_delta    = q0_for_delta;
    out.tail_frac       = tail_frac;
    out.tail_idx        = tail_idx;
    out.singular_values = s;
end