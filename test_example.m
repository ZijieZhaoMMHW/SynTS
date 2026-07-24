%% Example: generate synthetic SST with LIM-mode variability
%
% This script runs the example workflow included in the repository.
% Required files:
%   - gss_with_modes_fixed.m
%   - gss_old.m
%   - sst_here_example.mat
%   - lim_60k.mat

clear
clc

%% Settings

years_syn = 5000;        % number of synthetic years
seed = 10;               % random seed
acf_lags = 60;           % lags used for ACF diagnostics

obs_start = datenum(1982, 1, 1);
obs_end   = datenum(2022, 12, 31);
tnum_obs  = obs_start:obs_end;

%% Load data

load('sst_here_example.mat');     % provides sst_here
load('lim_60k.mat');              % provides lim_obs and lim_syn

dates_obs = datevec(tnum_obs);
is_leap_day = dates_obs(:,2) == 2 & dates_obs(:,3) == 29;
dates_obs_365 = dates_obs(~is_leap_day, :);

%% Convert monthly LIM modes to daily values

x_mode_obs = monthly_to_daily(lim_obs, dates_obs_365);

dates_syn = datevec(datenum(1,1,1):datenum(years_syn,12,31));
dates_syn = dates_syn(~(dates_syn(:,2) == 2 & dates_syn(:,3) == 29), :);

lim_syn = lim_syn(1:years_syn*12);
x_mode_syn = monthly_to_daily(lim_syn, dates_syn);

%% Generate synthetic SST

n_grid = size(sst_here, 1);
sst_syn = nan(n_grid, years_syn*365);
seas_mean_plot = [];
trend_obs_plot = [];

for i = 1:n_grid
    fprintf('Generating synthetic SST for grid point %d/%d...\n', i, n_grid);

    [sst_syn(i,:), ~, seas_mean, ~, trend_obs] = gss_with_modes_fixed( ...
        sst_here(i,:)', ...
        tnum_obs, ...
        x_mode_obs, ...
        x_mode_syn, ...
        seed + i - 1, ...
        years_syn);

    if i == 1
        seas_mean_plot = seas_mean;
        trend_obs_plot = trend_obs;
    end
end

%% Simple diagnostic plot for the first grid point

i_plot = 1;

sst_obs_365 = sst_here(i_plot, ~is_leap_day)';
obs_anom = sst_obs_365 - seas_mean_plot - trend_obs_plot;

syn_anom = remove_daily_climatology(sst_syn(i_plot,:)', years_syn);

acf_obs = autocorr(detrend(obs_anom), 'NumLags', acf_lags);
acf_syn = autocorr(detrend(syn_anom), 'NumLags', acf_lags);

figure('Color', 'w');
plot(0:acf_lags, acf_obs, 'k-', 'LineWidth', 2.5);
hold on
plot(0:acf_lags, acf_syn, 'r-', 'LineWidth', 2);
grid on
box on
xlabel('Lag (days)');
ylabel('Autocorrelation');
legend('Observed', 'Synthetic', 'Location', 'northeast');
title('Observed and synthetic SST anomaly memory');

%% Helper functions

function daily = monthly_to_daily(monthly, dates_365)
    month_days = [31 28 31 30 31 30 31 31 30 31 30 31];
    daily = nan(size(dates_365, 1), 1);

    for m = 1:12
        idx = dates_365(:,2) == m;
        vals = monthly(m:12:end);
        vals_daily = repmat(vals(:)', month_days(m), 1);
        daily(idx) = vals_daily(:);
    end
end

function anom = remove_daily_climatology(sst, years_syn)
    doy = repmat((1:365)', years_syn, 1);
    clim = nan(365, 1);

    for d = 1:365
        clim(d) = nanmean(sst(doy == d));
    end

    anom = sst - clim(doy);
end
