# Observation-constrained synthetic SST generator

This repository contains MATLAB code for generating long synthetic sea-surface
temperature (SST) records that preserve key statistical properties of the
observed ocean. The code was developed for a study on statistically defensible
detection of marine heatwaves and other threshold-based temperature extremes.

The central idea is to extend short observational records into long synthetic
records that can be used as a reference state for testing extreme-event
detection methods. The generated SST series preserve:

- the observed seasonal mean cycle,
- the observed seasonal cycle of variance,
- short- to subseasonal SST memory,
- large-scale climate-mode variability represented by LIM-derived modes,
- and a residual stochastic component fitted to the observed local SST record.

The code is intended as a lightweight research implementation rather than a
fully packaged toolbox.

## Main files

```text
gss_with_modes_fixed.m   Main function combining observed seasonal structure,
                         LIM-derived modes, and a stochastic residual.

gss_old.m                Residual SST generator used internally by
                         gss_with_modes_fixed.m.

test_example.m           Rough example showing how observed and synthetic LIM
                         modes can be expanded to daily resolution and passed
                         into the main function.
```

The example script expects local `.mat` files such as `sst_here_example.mat`
and `lim_60k.mat`. These data files are not part of the core code and should be
replaced with your own SST and mode data.

## Method overview

For each observed SST time series, `gss_with_modes_fixed.m` performs the
following steps:

1. Removes leap days and works on a 365-day calendar.
2. Fits the seasonal mean using harmonic regression.
3. Removes the seasonal mean and linear trend to obtain SST anomalies.
4. Estimates the day-of-year-dependent variance.
5. Splits the observed anomaly into a prescribed large-scale/modal component
   and a local residual.
6. Generates a long synthetic residual using `gss_old.m`.
7. Combines the synthetic residual with prescribed synthetic LIM-mode
   anomalies.
8. Rescales the combined anomaly so that its variance and 1-60 day
   autocorrelation are close to those of the observed anomaly.
9. Adds the observed seasonal mean back to obtain synthetic total SST.

In the paper, the synthetic SST records are used to provide a long reference
state against which marine-heatwave detection methods and baseline parameters
can be evaluated objectively.

## Requirements

The code is written for MATLAB.

Required MATLAB functionality:

- `datetime`, `datevec`, `datenum`
- `arima`, `estimate`, and `autocorr`
- `nanstd`, `nanmean`

These functions are commonly available with the Statistics and Machine
Learning Toolbox and Econometrics Toolbox. If `arima` or `autocorr` is missing,
check that the relevant MATLAB toolbox is installed.

## Function interface

```matlab
[sst_syn_total, sst_syn_res, seas_mean_obs, sigma_d_obs, trends_here] = ...
    gss_with_modes_fixed(sst_raw, tnum, x_mode_obs, x_mode_syn, SEED, YEARS_SYN)
```

### Inputs

```text
sst_raw
    Observed SST time series. Use a column vector or a vector that can be
    converted to a column vector. Units are usually deg C.

tnum
    MATLAB datenum vector corresponding to sst_raw. It must have the same
    length as sst_raw.

x_mode_obs
    Observed large-scale/modal SST anomaly used as the resolved climate-mode
    component. This should already be deseasonalized, detrended, and aligned
    to a 365-day calendar with leap days removed.

x_mode_syn
    Synthetic large-scale/modal SST anomaly, typically generated from a LIM or
    another climate-mode model. Its length must be YEARS_SYN * 365.

SEED
    Random seed used by the residual generator.

YEARS_SYN
    Number of synthetic years to generate. The output length is
    YEARS_SYN * 365.
```

### Outputs

```text
sst_syn_total
    Final synthetic total SST, including the seasonal mean and synthetic
    anomaly.

sst_syn_res
    Synthetic SST generated for the residual component only.

seas_mean_obs
    Fitted observed seasonal mean over the leap-day-filtered observed record.

sigma_d_obs
    Fitted day-of-year standard deviation of observed anomalies, length 365.

trends_here
    Fitted linear trend component removed from the observed deseasonalized SST.
```

## Important calendar convention

The code removes February 29 and assumes all synthetic years have 365 days.

This means:

- `sst_raw` and `tnum` may include leap days; they are filtered internally.
- `x_mode_obs` should already be on the leap-day-removed 365-day calendar.
- `x_mode_syn` must have length `YEARS_SYN * 365`.

If your modal input still contains leap days, either remove them before calling
the function or adapt the commented filtering line in `gss_with_modes_fixed.m`.

## Minimal usage pattern

```matlab
% Observed SST and corresponding dates
load('sst_here_example.mat');     % should provide sst_here or similar
tnum_obs = datenum(1982,1,1):datenum(2022,12,31);

% Observed and synthetic modal anomalies
% These should be prepared outside the function, for example from LIM output.
load('lim_60k.mat');              % should provide lim_obs and lim_syn or similar

YEARS_SYN = 5000;
SEED = 10;

% Example placeholders. Replace these with your actual mode-processing code.
x_mode_obs = sstlim_obs(:);       % length should match observed non-leap days
x_mode_syn = sstlim_daily(:);     % length should be YEARS_SYN * 365

[sst_syn_total, sst_syn_res, seas_mean_obs, sigma_d_obs, trends_here] = ...
    gss_with_modes_fixed(sst_here(:), tnum_obs, x_mode_obs, x_mode_syn, SEED, YEARS_SYN);
```

For a rough end-to-end example, see `test_example.m`.

## Preparing daily LIM-mode inputs

The example script assumes monthly LIM-mode anomalies and expands them to daily
resolution by repeating each monthly value for the number of days in that
month. In outline:

```matlab
mnum = [31 28 31 30 31 30 31 31 30 31 30 31];

for month = 1:12
    idx_here = date_used(:,2) == month;
    ts_here = repmat(lim_monthly(month:12:end), mnum(month), 1);
    lim_daily(idx_here) = ts_here(:);
end
```

You can replace this with a more sophisticated interpolation or mode-generation
approach if needed. The key requirement is that the observed and synthetic mode
vectors are aligned to the same 365-day calendar convention used by the SST
generator.

## Quick diagnostic check

After generating synthetic SST, compare the observed and synthetic anomaly
autocorrelation functions:

```matlab
acf_obs = autocorr(detrend(obs_anom(:)), 'NumLags', 60);
acf_syn = autocorr(detrend(syn_anom(:)), 'NumLags', 60);

figure
plot(acf_obs, 'LineWidth', 2); hold on
plot(acf_syn, 'LineWidth', 2);
legend('Observed', 'Synthetic');
xlabel('Lag (days)');
ylabel('Autocorrelation');
```

The rough example in `test_example.m` includes a similar check comparing
observed and synthetic ACFs.

## Common pitfalls

### Dimension mismatch between SST and mode inputs

If MATLAB reports incompatible array sizes near:

```matlab
x_res_obs = x_total - x_mode_obs;
```

then `x_mode_obs` is probably not aligned to the leap-day-removed observed SST
record. Remove February 29 from the modal input before calling the function.

### Synthetic mode length mismatch

`x_mode_syn` must have exactly `YEARS_SYN * 365` elements. If it is generated
from monthly modes, expand it to daily resolution first.

### Missing `gss_old.m`

`gss_with_modes_fixed.m` calls `gss_old.m` internally. Keep both files in the
same MATLAB path.

### Very long synthetic records

Generating thousands of synthetic years can be memory-intensive, especially
when looping over many grid cells. For testing, start with a smaller value such
as:

```matlab
YEARS_SYN = 100;
```

Then increase to the desired production length once the workflow is validated.

## Suggested repository structure

```text
.
├── README.md
├── gss_with_modes_fixed.m
├── gss_old.m
├── test_example.m
├── data/
│   └── README.md              # describe where example data can be obtained
└── examples/
    └── test_example.m
```

Large `.mat` files should generally not be committed directly to GitHub. If
example data are needed, consider hosting them separately and linking to them
from the repository.

## Citation

If you use this code, please cite the associated manuscript:

```text
Zhao et al. Detecting Temperature Extremes in a Statistically Optimal Way.
Manuscript in preparation/submitted.
```

Update this citation once the final publication information is available.
