function [sst_syn,seas_mean,syn_var,sigma_d] = gss_old(sst_raw,tnum,SEED,YEARS_SYN)

% ---------------------- Config ----------------------
P_ANNUAL = 365.2425;
H_MEAN = 3;      % # of harmonics for seasonal mean
H_LOGVAR = 3;    % # of harmonics for log-variance
%YEARS_SYN = 5000;
DAYS_PER_YR = 365;    % drop Feb 29; assume 365-day years
N_SYN = YEARS_SYN * DAYS_PER_YR;
%SEED = 20251029;

rng(SEED)
% sst_raw = double(squeeze(ncread('sst_daily_1D.nc','sst')));
% sst_raw = sst_raw(1:15706);
% tnum = datenum(1982,1,1):datenum(2024,12,31);

dates = datevec(tnum);
dates(:,1)=2001;
keep = ~(dates(:,2) == 2 & dates(:,3) == 29);
sst = sst_raw(keep);
dates = dates(keep,:);
t = (0:length(sst)-1)';
doy=day(datetime(dates),'dayofyear');

Xs = design_harmonics(t, P_ANNUAL, H_MEAN);
beta = Xs \ sst;
seas_mean = Xs * beta;

x0 = sst - seas_mean;
tc = t - mean(t);
Xt = [ones(size(tc)), tc];
btr = Xt \ x0;
x = x0 - (Xt * btr);

% 3) Seasonal variance via log-variance regression (2 harmonics)
raw_var = nan(DAYS_PER_YR, 1);
for d = 1:DAYS_PER_YR
    sel = x(doy == d);
    if length(sel) >= 2
        raw_var(d) = var(sel, 0);
    end
end
gvar = var(x, 0);
raw_var(~isfinite(raw_var)) = gvar;

dgrid = (1:DAYS_PER_YR)';
Xv = design_harmonics(dgrid, P_ANNUAL, H_LOGVAR);
syn_var=Xv * (Xv \log(raw_var));
gamma = Xv \ log(raw_var);
sigma_d = exp(0.5 * (Xv * gamma));        % Ïƒ(d) for DOY 1..365
sigma_t = sigma_d(doy);                   % map to each observed day

% 4) Fit ARMA(1,1) to Ïƒ-standardised anomalies
z = x ./ sigma_t;

% Fit ARMA(1,1) using MATLAB's arima
%fprintf('Fitting ARMA(1,1) model...\n');
Mdl = arima('ARLags', 1, 'MALags', 1, 'Constant', 0);
EstMdl = estimate(Mdl, z, 'Display', 'off');
phi = EstMdl.AR{1};
theta = EstMdl.MA{1};
sig2 = EstMdl.Variance;
std_z = std(z, 0);

% 5) Optimise slow memory (rho, f) vs observed ACF(10-60 d)
N_ACF = 60;
acf_obs = autocorr(x, 'NumLags', N_ACF);
lags = 0:N_ACF;
mask = (lags >= 10) & (lags <= 60);

% Generate fast component for optimization
y_fast = simulate_arma11(phi, theta, length(z), sqrt(sig2));
y_fast = y_fast / std(y_fast, 0);

% Objective function
obj_func = @(params) objective(params(1), params(2), y_fast, std_z, ...
    sigma_t, mask, acf_obs, N_ACF);

% Grid search
best_val = inf;
best_rho = 0;
best_f = 0;

% Coarse search
for rho = 0.90:0.003:0.999
    for f = 0.05:0.05:0.85
        val = obj_func([rho, f]);
        if val < best_val
            best_val = val;
            best_rho = rho;
            best_f = f;
        end
    end
end

% Fine search around best values
rho_range = max(0.90, best_rho-0.003):0.001:min(0.9999, best_rho+0.003);
f_range = max(0.01, best_f-0.08):0.02:min(0.99, best_f+0.08);

for rho = rho_range
    for f = f_range
        val = obj_func([rho, f]);
        if val < best_val
            best_val = val;
            best_rho = rho;
            best_f = f;
        end
    end
end

rho_opt = best_rho;
f_opt = best_f;

% 6) Long simulation: 50,000 years of SST

yA = simulate_arma11(phi, theta, N_SYN, sqrt(sig2));
yA = yA / std(yA, 0);

yS = simulate_ar1(rho_opt, N_SYN);
yS = yS / std(yS, 0);

y_mix = sqrt(1 - f_opt) * yA + sqrt(f_opt) * yS;
y_mix = y_mix * std_z;  % match z scale

sigma_long = repmat(sigma_d, YEARS_SYN, 1);
x_syn = y_mix .* sigma_long;  % synthetic anomalies in Â°C

harm_1yr = design_harmonics((0:DAYS_PER_YR-1)', P_ANNUAL, H_MEAN) * beta;
sst_syn = repmat(harm_1yr, YEARS_SYN, 1) + x_syn;
doy_long = mod(0:N_SYN-1, DAYS_PER_YR)' + 1;
end

function X = design_harmonics(t, period, H)
    % Return [1, cos(2Ï€kt/P), sin(2Ï€kt/P)] for k=1..H
    t = double(t(:));
    w = 2 * pi / period;
    X = ones(length(t), 1 + 2*H);
    
    for k = 1:H
        X(:, 2*k) = cos(k * w * t);
        X(:, 2*k+1) = sin(k * w * t);
    end
end

function y = simulate_arma11(phi, theta, n, noise_std)
    % ARMA(1,1): y_t = phi*y_{t-1} + e_t + theta*e_{t-1}
    e = randn(n, 1) * noise_std;
    y = zeros(n, 1);
    y(1) = e(1);
    
    for t = 2:n
        y(t) = phi * y(t-1) + e(t) + theta * e(t-1);
    end
end

function s = simulate_ar1(rho, n)
    % AR(1) with stationary initial variance
    eta = randn(n, 1);
    s = zeros(n, 1);
    s(1) = eta(1) / sqrt(max(1e-12, 1 - rho^2));
    
    for t = 2:n
        s(t) = rho * s(t-1) + eta(t);
    end
end

function val = objective(rho, f, y_fast, std_z, sigma_t, mask, acf_obs, N_ACF)
    % Objective function for parameter optimization
    s_short = simulate_ar1(rho, length(y_fast));
    s_unit = s_short / std(s_short, 0);
    y = sqrt(1 - f) * y_fast + sqrt(f) * s_unit;
    y = y * std_z .* sigma_t;  % back to Â°C anomalies for ACF
    ac = autocorr(y, 'NumLags', N_ACF);
    val = sum((ac(mask) - acf_obs(mask)).^2);
end