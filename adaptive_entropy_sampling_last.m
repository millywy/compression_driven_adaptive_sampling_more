% adaptive_entropy_sampling.m
% Entropy-driven adaptive sampling atop the WFPV (TBME2017) pipeline.
% Two modes: fs_hi=25 Hz, fs_lo=6.25 Hz. Raw data assumed at fs0=125 Hz.

clearvars; % Clears all variables from the workspace
%clc;       % Clears the command window
close all; % Closes all figure windows

%% Configuration
fs0 = 125;
fs_hi = 25;
fs_lo = 50;
fs_proc = 125;             % internal WFPV processing rate (matches original)
FORCE_LOW = true;         % set true to lock controller to LOW for baseline comparison
fs_acc = 25;              % fixed-rate control stream for ACC
FFTres = 1024;
WFlength = 15;            % Wiener averaging length (frames)
CutoffFreqHzBP = [0.4 4];     % bandpass at 125 Hz before decimation (Nyquist at 3.125 Hz)
CutoffFreqHzSearch = [1 3];   % HR search band (Hz)
window_sec = 8;
step_sec = 2;

% Adaptive entropy thresholds / hysteresis (tunable, should be parametrizable)
nbits_entropy = 2;            % quantization for entropy proxy
hi_hold = 20;                 % min windows to stay high after going HIGH
Th_hi = 0.15;                  % unstable = above this
N_look_back = 7;               % windows to look back for stable entropy
N_unstable = 2;                % unstable windows to go HIGH
Th_low = 0.13;                 % stable = below this
N_stable  = 5;                % min stable windows in high before going LOW


IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

%% Metrics containers
logRec = struct([]);
myError = nan(1, numel(IDData));

%% Process each recording
for idnb = 1:numel(IDData)
    % Load data and pick PPG/ACC channels (dataset-specific)
    load(['Data/' IDData{idnb}], 'sig');
    if idnb > 13
        ch = [1 2 3 4 5];
    else
        ch = [2 3 4 5 6];
    end

    sig_raw = sig(ch, :);

    % Sliding window setup (seconds -> samples)
    window = window_sec * fs0;
    step   = step_sec   * fs0;
    windowNb = floor((size(sig_raw,2)-window)/step) + 1;

    if windowNb < 1
        % Skip very short recordings
        BPM_est = [];
        Hacc = [];
        ACCmag25 = [];
        return;
    end

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
    fs_cur = fs_lo;                 % start low
    fs_next = fs_lo;
    state_mode = "LOW";             % start low
    hi_timer = 0;                   % tracks how long we must stay HIGH

    
    for i = 1:windowNb
        % Extract current frame from raw data
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        curDataRaw = sig_raw(:, curSegment);

        % Run one-frame HR est with internal resampling (match baseline helper);
        % 1)load in correct mode state
        if fs_cur == fs_hi
            state_in = state_hi;
        else
            state_in = state_lo;
        end

        % 2) process HR in this frame using current mode
        [BPM_est(i), state_out] = WFPV_bb.simple_hr_one_frame_last(curDataRaw, fs0, fs_cur, fs_proc, FFTres, WFlength, CutoffFreqHzSearch, state_in, i, BPM_est, idnb, CutoffFreqHzBP);
        % [BPM_est(i), state_out] = WFPV_bb.wfpv_one_frame_last(curDataRaw, fs0, fs_cur, fs_proc, FFTres, WFlength, CutoffFreqHzSearch, state_in, i, BPM_est, idnb, CutoffFreqHzBP);

        % 3) Save back mode state
        if fs_cur == fs_hi
            state_hi = state_out;
        else
            state_lo = state_out;
        end

        %% Filter at original rate (reuse per recording)
        [b125, a125] = butter(4, [0.3 8]/(fs0/2), 'bandpass');
        % 4) Filter ACC at 125 Hz before downsampling
        curDataFilt2 = curDataRaw;
        for c = 3:size(curDataRaw,1)
            curDataFilt2(c,:) = filter(b125, a125, curDataRaw(c,:));
        end

        % 5) ACC control stream at fixed 25 Hz for entropy calculation
        curAcc_resampled = do_resample_last(curDataFilt2, fs0, fs_acc);
        ACCmag25 = sqrt(curAcc_resampled(1,:).^2 + curAcc_resampled(2,:).^2 + curAcc_resampled(3,:).^2);
        ACCmag25_log(i) = mean(ACCmag25);
        Hacc(i) = entropy_proxy_hist(ACCmag25, nbits_entropy);
        if i == 1
            dHacc(i) = 0;
        else
            dHacc(i) = Hacc(i) - Hacc(i-1);
        end

        % % 5) ACC control stream at fixed 125 Hz for entropy calculation
        % ACCmag25 = sqrt(curDataFilt2(1,:).^2 + curDataFilt2(2,:).^2 + curDataFilt2(3,:).^2);
        % ACCmag25_log(i) = mean(ACCmag25);
        % Hacc(i) = entropy_proxy_hist(ACCmag25, nbits_entropy);
        % if i == 1
        %     dHacc(i) = 0;
        % else
        %     dHacc(i) = Hacc(i) - Hacc(i-1);
        % end

        % 6) Adaptive sampling controller logic (simple hysteresis on entropy jumps)
        i0 = max(1, i - N_look_back + 1);
        recent = dHacc(i0:i);
        enough_hist = (i >= N_look_back);

        % State machine controller
        if FORCE_LOW
            state_mode = "LOW"; fs_next = fs_lo;
        else
            switch state_mode
                case "LOW"
                    fs_cur = fs_lo;
                    if enough_hist
                        n_bad = sum(abs(recent) > Th_hi);         % count unstable windows in look-back
                        if n_bad >= N_unstable
                            state_mode = "HIGH";
                            hi_timer = hi_hold;              % force HIGH for at least hi_hold windows
                            fs_next = fs_hi;
                        end
                    end
                case "HIGH"
                    fs_cur = fs_hi;
                    % Decrement hold timer while in HIGH
                    if hi_timer > 0
                        hi_timer = hi_timer - 1;
                    end

                    % After hold expires, exit only if ALL last 5 are stable (< Th_low)
                    if enough_hist && hi_timer == 0
                        n_good = sum(abs(recent) < Th_low);   % count stable in last N_look_back
                        if n_good >= N_stable
                            state_mode = "LOW";
                            fs_next = fs_lo;
                        end
                    end
                
            end
        end
        FsUsed(i) = fs_next;
        fs_cur = fs_next;
    end
    

    % Ground truth and error
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


    % Entropy summary stats
    stats.Hacc.median = median(Hacc);
    stats.Hacc.p25 = prctile(Hacc,25);
    stats.Hacc.p75 = prctile(Hacc,75);
    stats.Hacc.max = max(Hacc);
    logRec(idnb).entropyStats = stats;
end

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
state.FreqRange = [];
end
