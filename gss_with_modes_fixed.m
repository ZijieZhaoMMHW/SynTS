function [sst_syn_total, sst_syn_res, seas_mean_obs, sigma_d_obs,trends_here] = ...
    gss_with_modes_fixed(sst_raw, tnum, x_mode_obs, x_mode_syn, SEED, YEARS_SYN)
% Revised version: generate synthetic SST series with the observed seasonal cycle/variance/memory + LIM modes
%
% Key fixes:
%   1. Unified standardization
%   2. Correct variance allocation
%   3. Consistent seasonal variance handling
%
% Inputs:
%   sst_raw    : Observed SST (column vector, length N_obs)
%   tnum       : Corresponding datenum values (same length)
%   x_mode_obs : Observed modal anomalies (°C, with seasonal mean and trend removed)
%   x_mode_syn : Synthetic modal anomalies (°C, length N_syn = YEARS_SYN*365)
%   SEED       : Random seed
%   YEARS_SYN  : Number of synthetic years
%
% Outputs:
%   sst_syn_total : Final synthetic total SST (mode + residual)
%   sst_syn_res   : Synthetic SST for the residual component only
%   seas_mean_obs : Observed seasonal mean
%   sigma_d_obs   : Observed DOY standard deviation (365x1)
%   diagnostics   : Diagnostic information structure

    % ============ Configuration ============
    P_ANNUAL   = 365.2425;
    H_MEAN     = 3;
    H_LOGVAR   = 3;
    DAYS_PER_YR = 365;
    N_SYN      = YEARS_SYN * DAYS_PER_YR;

    % ============ Step 1: Preprocess observational data ============
    
    dates = datevec(tnum);
    dates(:,1) = 2001;
    keep = ~(dates(:,2) == 2 & dates(:,3) == 29);
    sst = sst_raw(keep);
    %x_mode_obs = x_mode_obs(keep);  % Apply the same filtering
    dates = dates(keep,:);
    t = (0:length(sst)-1)';
    doy = day(datetime(dates),'dayofyear');

    % 1.1 Seasonal mean
    Xs = design_harmonics_local(t, P_ANNUAL, H_MEAN);
    beta = Xs \ sst;
    seas_mean_obs = Xs * beta;

    % 1.2 Remove seasonality and trend
    x0 = sst - seas_mean_obs;
    tc = t - mean(t);
    Xt = [ones(size(tc)), tc];
    btr = Xt \ x0;
    trends_here=(Xt * btr);
    x_total = x0 - (Xt * btr);  % Total anomaly (°C)

    % 1.3 Estimate DOY variance sigma(d)
    raw_var = nan(DAYS_PER_YR,1);
    for d = 1:DAYS_PER_YR
        sel = x_total(doy == d);
        if numel(sel) >= 2
            raw_var(d) = var(sel,0);
        end
    end
    gvar = var(x_total,0);
    raw_var(~isfinite(raw_var)) = gvar;

    dgrid = (1:DAYS_PER_YR)';
    Xv = design_harmonics_local(dgrid, P_ANNUAL, H_LOGVAR);
    gamma = Xv \ log(raw_var);
    logvar_fit = Xv * gamma;
    sigma_d_obs = exp(0.5 * logvar_fit);
    sigma_t_obs = sigma_d_obs(doy);
    z = x_total ./ sigma_t_obs;
    std_z = nanstd(z);
    

    % ============ Step 2: Decompose into mode + residual ============
    
    % Ensure x_mode_obs is already deseasonalized and detrended
    % If not, preprocess it first
    x_mode_obs = x_mode_obs(:);
    
    % Compute the residual
    x_res_obs = x_total - x_mode_obs;

    % ============ Step 3: Generate the residual using GSS ============
    
    % Construct the observed series for the residual-only component
    sst_res_raw = seas_mean_obs + (Xt * btr) + x_res_obs;

    % Call gss_old
    [sst_syn_res, ~, ~, sigma_d_res] = ...
        gss_old(sst_res_raw, tnum(keep), SEED, YEARS_SYN);

    % ============ Step 4: Extract residual anomalies ============
    
    doy_long = repmat((1:DAYS_PER_YR)', YEARS_SYN, 1);
    
    % Compute the seasonal mean of the synthetic residual
    seas_res_syn = nan(DAYS_PER_YR,1);
    for d = 1:DAYS_PER_YR
        seas_res_syn(d) = mean(sst_syn_res(doy_long == d));
    end
    seas_res_long = seas_res_syn(doy_long);
    
    % Residual anomaly
    x_res_syn = sst_syn_res - seas_res_long;

    % ============ Step 5: Trial mixing ============
%     std_mode = nanstd(x_mode_syn(:));
%     std_res = nanstd(x_res_syn(:));
    
%     s_mode_syn=x_mode_syn(:)./std_mode;
%     s_res_syn=x_res_syn(:)./std_res;
 
    s_mode_syn=x_mode_syn(:);
    s_res_syn=x_res_syn(:);
    
    for i=1:365
        idx_here=doy_long==i;
        s_mode_syn(idx_here)=x_mode_syn(idx_here)./nanstd(x_mode_syn(idx_here));
        s_res_syn(idx_here)=x_res_syn(idx_here)./nanstd(x_res_syn(idx_here));
    end
    
    acf_obs=autocorr(x_total,'NumLags',60);
    
    f_used=0:0.03:1;
    for f=1:length(f_used)
        s_mix=sqrt(1-f_used(f))*s_mode_syn(1:length(z))+sqrt(f_used(f))*s_res_syn(1:length(z));
        s_mix=s_mix./nanstd(s_mix).*std_z.*sigma_t_obs;
        ac = autocorr(s_mix, 'NumLags', 60);
        val(f) = sum((ac(1:60) - acf_obs(1:60)).^2);
    end
    [~,loc]=nanmin(val);
    s_mix=sqrt(1-f_used(loc))*s_mode_syn(:)+sqrt(f_used(loc))*s_res_syn(:);
    x_total_syn=s_mix./nanstd(s_mix).*std_z.*repmat(sigma_d_obs,YEARS_SYN,1);

    % ============ Step 6: Recombine ============
    
    % Synthetic total anomaly
    %x_total_syn = x_mode_syn(:) + x_res_syn(:);
    
    sst_syn_total = repmat((seas_mean_obs(1:365)),YEARS_SYN,1) + x_total_syn;

    
    
end


% ---- Helper function ----
function X = design_harmonics_local(t, period, H)
    t = double(t(:));
    w = 2*pi / period;
    X = ones(length(t), 1 + 2*H);
    for k = 1:H
        X(:, 2*k)   = cos(k*w*t);
        X(:, 2*k+1) = sin(k*w*t);
    end
end
