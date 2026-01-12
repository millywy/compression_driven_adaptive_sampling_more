## added to overall configuration
%% Configuration
fs0 = 125;               % original sampling rate
fs_hi = 25;
fs_lo = 12.5;
fs_proc = 25;             % internal WFPV processing rate (matches original)
FORCE_LOW = true;         % set true to lock controller to LOW for baseline comparison
fs_acc = 25;              % fixed-rate control stream for ACC
<!-- CutoffFreqHzBP = [0.4 4];     % bandpass at 125 Hz before decimation (Nyquist at 3.125 Hz)
CutoffFreqHzSearch = [1 3];   % HR search band (Hz) -->

% Adaptive entropy thresholds / hysteresis (tunable, should be parametrizable)
nbits_entropy = 2;            % quantization for entropy proxy
hi_hold = 20;                 % min windows to stay high after going HIGH
Th_hi = 0.15;                  % unstable = above this
N_look_back = 7;               % windows to look back for stable entropy
N_unstable = 2;                % unstable windows to go HIGH
Th_low = 0.13;                 % stable = below this
N_stable  = 5;                % min stable windows in high before going LOW

%% Metrics containers
logRec = struct([]);
myError = nan(1, numel(IDData));

## added to per recording codes
BPM_est = zeros(1, windowNb);
FsUsed  = zeros(1, windowNb);
Hacc    = zeros(1, windowNb);          % ACC entropy at fixed 25 Hz
dHacc = zeros(1, windowNb);            % ACC entropy diff (jump metric)
ACCmag25_log = zeros(1, windowNb);     % mean ACC magnitude per window

% Mode states (keep WF history per mode so HIGH/LOW do not mix)
state_hi = init_mode_state();
state_lo = init_mode_state();
state_in = init_mode_state();
state_out = init_mode_state();
fs_adc = fs_lo;                 % start low
fs_next = fs_lo;
state_mode = "LOW";             % start low
hi_timer = 0;                   % tracks how long we must stay HIGH


## added to per frame codes
% Run one-frame HR est with internal resampling (match baseline helper);

% 1) bandpass at 125 Hz, downsample to fs_adc, then resample to fs_proc
[b,a] = butter(4, bpHz/(fs0/2), 'bandpass');
curDataFilt = zeros(size(curDataRaw));
for c = 1:size(curDataRaw,1)
    curDataFilt(c,:) = filter(b,a,curDataRaw(c,:));
end

% downsample to fs_adc if needed
if abs(fs0 - fs_adc) < eps
    curData_adc = curDataFilt; fs = fs_proc;
else
    curData_adc = do_resample_last(curDataFilt, fs0, fs_adc);
    fs = fs_proc;
end

% resample to internal if needed
if abs(fs_adc - fs_proc) < eps
    curData = curData_adc; fs = fs_proc;
else
    curData = do_resample_last(curData_adc, fs_adc, fs_proc);
    fs = fs_proc;
end

% 2)load in correct mode state
if fs_adc == fs_hi
    state_in = state_hi;
else
    state_in = state_lo;
end

W1_FFTi = state_in.W1_FFTi;
W11_FFTi = state_in.W11_FFTi;
W2_FFTi = state_in.W2_FFTi;
W21_FFTi = state_in.W21_FFTi;
PPG_ave_FFTpr = state_in.prevFFT;
rangeIdx = state_in.rangeIdx;

% 3) process HR in this frame using current mode, curData
<!-- [BPM_est(i), state_out] = WFPV_bb.simple_pd_hr_one_frame(curData); -->
take out their filtering and downsampling
edit places where pr mode saving is needed- eg in vocoder, arrays not same lenght bcs 6.25Hz FFT has ~150bins in search range vs 25Hz FFT as ~600bins in search range

% 3) Save back mode state

state_out.W1_FFTi = W1_FFTi;
state_out.W11_FFTi = W11_FFTi;
state_out.W2_FFTi = W2_FFTi;
state_out.W21_FFTi = W21_FFTi;
state_out.prevFFT = PPG_ave_FFTpr;
state_out.rangeIdx = rangeIdx;

if fs_adc == fs_hi
    state_hi = state_out;
else
    state_lo = state_out;
end

% 3.5) compute entropy

% 4) frequency decision 2 state FSM
[fs_next, next_state_mode] = fadc_decision(curDataRaw, ACCmag25_log, Hacc, dHacc, FsUsed, fs_adc, state_mode);
fs_adc = fs_next;
state_mode = next_state_mode;


% 5)Ground truth and error
if idnb > 13
    load(['Data/True' IDData{idnb}(5:end)], 'BPM0');
else
    load(['Data/' IDData{idnb} '_BPMtrace'], 'BPM0');
end
frames = min(length(BPM_est), length(BPM0));
myError(idnb) = mean(abs(BPM0(1:frames) - BPM_est(1:frames)'));

% Store log
logRec(idnb).ID = IDData{idnb};
logRec(idnb).FsUsed = FsUsed;
logRec(idnb).Hacc = Hacc;
logRec(idnb).dHacc = dHacc;
logRec(idnb).BPM_est = BPM_est(1:frames);
logRec(idnb).BPM0 = BPM0(1:frames);
logRec(idnb).ACC_mag_25 = ACCmag25_log(1:frames);


## added to the end
%% Aggregate metrics
MAE_all = mean(myError, 'omitnan');
MAE_train = mean(myError(1:12), 'omitnan');
MAE_test = mean(myError(13:end), 'omitnan');
fprintf('\n=== Adaptive Entropy Sampling Results ===\n');
fprintf('Err12=%2.2f, Err11=%2.2f, ErrAll=%2.2f\n', MAE_train, MAE_test, MAE_all);
fprintf('Individual recording errors (BPM):\n');
fprintf(' ');
fprintf('%4.2f ', myError);
fprintf('\n');

% Bland-Altman and correlation (aggregate all frames)
fullBPM0 = [];
fullBPM = [];
for rr = 1:numel(logRec)
    fullBPM0 = [fullBPM0, logRec(rr).BPM0(:)'];
    fullBPM  = [fullBPM,  logRec(rr).BPM_est(:)'];
end
fprintf('Generating Bland-Altman plot...\n');
[~, figBA] = BlandAltman(fullBPM0', fullBPM', {'Ground truth HR','Estimated HR'});
if exist('sgtitle','file') && ~isempty(figBA), figure(figBA); sgtitle('Bland-Altman (Adaptive Entropy Sampling)'); end
tmp = corrcoef(fullBPM0, fullBPM);
fprintf('Overall correlation coefficient: %.4f\n', tmp(1,2));

% Mode usage stats
all_fs_used = [logRec.FsUsed];
pct_hi = 100 * mean(all_fs_used==fs_hi);
pct_lo = 100 * mean(all_fs_used==fs_lo);
fprintf('Mode usage: high=%.1f%% low=%.1f%%\n', pct_hi, pct_lo);

% Persist summary logs for later analysis
save('adaptive_entropy_logs.mat', 'logRec', 'myError', 'MAE_all', 'MAE_train', 'MAE_test');

%% Selected recordings for comparison
selRecs = [10 17 19];
for idx = 1:numel(selRecs)
    r = selRecs(idx);
    if r <= numel(logRec) && ~isempty(logRec(r).BPM0)
        win_count = numel(logRec(r).Hacc);
        frames = 1:win_count;

        % HR comparison and controller signals
        figure;
        plot(logRec(r).BPM0,'ro'); hold on; 
        plot(logRec(r).BPM_est,'o','Color','blue'); 
        xlabel('Time (frames)'); ylabel('HR (BPM)'); legend({'Ground truth','Estimates'});
        yyaxis right;
        plot(logRec(r).Hacc, '--', 'Color', 'green', 'LineWidth', 1.0, 'DisplayName', 'Hacc (25 Hz ACC)');
        plot(logRec(r).dHacc, '-', 'Color', [0 1 1], 'LineWidth', 1.0, 'DisplayName', 'dHacc (25 Hz ACC)');
        yline(0.15, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'Ref 0.15');
        yline(-0.15, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'Ref -0.15');
        ylabel('Hacc / dHacc');
        yyaxis left;
        stairs(frames, logRec(r).FsUsed, '-'); ylabel('FsUsed (Hz)');
        title(sprintf('Recording %d (Adaptive Entropy Sampling)', r));
        grid on;

    end
end

%% Helpers
function state = init_mode_state()
% Container for per-mode FFT/Wiener history so HIGH/LOW stay independent
state.W1_FFTi = [];
state.W11_FFTi = [];
state.W2_FFTi = [];
state.W21_FFTi = [];
state.prevFFT = [];
state.rangeIdx = [];
end
