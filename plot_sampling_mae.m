% plot_sampling_mae.m
% Quick summary plot of HR MAE (bpm) versus effective PPG sampling rate (Hz)
% using the numbers you reported for uniform baselines and four adaptive
% controllers (HIST, ARITH-O1, LZ, ACC magnitude controller).
%
% Effective rate:
%   R_eff = p_high * f_high + (1 - p_high) * f_low
% Assumes f_high = 25 Hz and f_low = 6.25 Hz for adaptive runs.

clear;

f_high = 25;
f_low  = 6.25;

% Data you provided (MAE_all) -------------------------
% Uniform baselines (actual ADC rate)
uniform_fs  = [25    12.5   6.25];
uniform_mae = [3.01  6.02   6.50];

% Adaptive controllers (p_high in %, MAE_all)
adapt_labels = {'Entropy HIST','Entropy ARITH-O1','Entropy LZ','ACC mag controller'};
p_high_pct   = [55.0,         58.7,               73.4,        68.1];
mae_adapt    = [4.55,         4.32,               3.45,        4.74];

% Compute effective sampling for adaptive runs
p_high = p_high_pct / 100;
fs_eff_adapt = p_high * f_high + (1 - p_high) * f_low;

% Consolidate for plotting
labels = [ ...
    "Uniform 25 Hz", ...
    "Uniform 12.5 Hz", ...
    "Uniform 6.25 Hz", ...
    adapt_labels ...
    ];
fs_eff = [uniform_fs, fs_eff_adapt];
mae    = [uniform_mae, mae_adapt];

% Plot
figure;
plot(fs_eff, mae, 'o-', 'LineWidth', 1.4); grid on;
xlabel('Effective sampling frequency (Hz)');
ylabel('MAE (BPM)');
title('WFPV HR MAE vs Effective Sampling Rate');

% Add text labels near points
for k = 1:numel(fs_eff)
    text(fs_eff(k)+0.1, mae(k), labels(k), 'Interpreter','none');
end

% Print summary table
fprintf('\nSummary: MAE vs effective sampling rate\n');
fprintf('%-25s  fs_eff(Hz)   MAE(BPM)\n', 'Label');
for k = 1:numel(fs_eff)
    fprintf('%-25s  %8.2f    %6.2f\n', labels(k), fs_eff(k), mae(k));
end

% -------- Correlation coefficient vs effective sampling rate --------
% You reported corrcoef values for each run; add them here.
% For adaptive runs:
corr_adapt = [0.9135, 0.9226, 0.9543, 0.8930]; % HIST, ARITH-O1, LZ, ACC mag
% For uniform runs:
corr_uniform = [0.9507, 0.8145, 0.8099];       % 25, 12.5, 6.25 Hz

corr_vals = [corr_uniform, corr_adapt];

figure;
plot(fs_eff, corr_vals, 's-', 'LineWidth', 1.4); grid on;
xlabel('Effective sampling frequency (Hz)');
ylabel('Correlation coefficient');
title('Correlation vs Effective Sampling Rate');
for k = 1:numel(fs_eff)
    text(fs_eff(k)+0.1, corr_vals(k), labels(k), 'Interpreter','none');
end

fprintf('\nSummary: Corrcoef vs effective sampling rate\n');
fprintf('%-25s  fs_eff(Hz)   Corr\n', 'Label');
for k = 1:numel(fs_eff)
    fprintf('%-25s  %8.2f    %6.4f\n', labels(k), fs_eff(k), corr_vals(k));
end

% -------- Uniform-only overlay: MAE (left) and Corr (right) --------
figure;
yyaxis left;
plot(uniform_fs, uniform_mae, 'o-', 'LineWidth', 1.4, 'DisplayName','MAE');
ylabel('MAE (BPM)');

yyaxis right;
plot(uniform_fs, corr_uniform, 's--', 'LineWidth', 1.4, 'DisplayName','Corr');
ylabel('Correlation coefficient');

xlabel('Uniform sampling frequency (Hz)');
title('Baseline MAE and Corr under Uniform Decimation');
grid on;
legend('Location','best');
