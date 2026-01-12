function [fs_next, next_state_mode, hi_timer_out, hacc_val, dhacc_val, accmag_val] = ...
    fadc_decision(curDataRaw, fs0, fs_lo, fs_hi, state_mode, i, hi_timer_in, Hacc, dHacc)

% fs0 = 125;               % original sampling rate
% fs_hi = 25;
% fs_lo = 6.25;
% fs_proc = 25;             % internal WFPV processing rate (matches original)
FORCE_LOW = false;         % set true to lock controller to LOW for baseline comparison
% fs_next = fs_adc;        % default to current rate
% % fs_acc = 25;              % fixed-rate control stream for ACC
% % CutoffFreqHzBP = [0.4 4];     % bandpass at 125 Hz before decimation (Nyquist at 3.125 Hz)
% % CutoffFreqHzSearch = [1 3];   % HR search band (Hz) 


%Adaptive entropy thresholds / hysteresis (tunable, should be parametrizable)
nbits_entropy = 2;            % quantization for entropy proxy
hi_hold = 20;                 % min windows to stay high after going HIGH
Th_hi = 0.15;                  % unstable = above this
N_look_back = 7;               % windows to look back for stable entropy
N_unstable = 2;                % unstable windows to go HIGH
Th_low = 0.13;                 % stable = below this
N_stable  = 5;                % min stable windows in high before going LOW

    %% Process ACC for entropy calculation

% Filter ACC at 125 Hz before entropy calculation
[b125, a125] = butter(4, [0.3 8]/(fs0/2), 'bandpass'); %can modify freq range here
curDataFilt2 = curDataRaw;
for c = 3:size(curDataRaw,1)
    curDataFilt2(c,:) = filter(b125, a125, curDataRaw(c,:));
end

% % ACC control stream at fixed 25 Hz for entropy calculation
% curAcc_resampled = do_resample_last(curDataFilt2, fs0, fs_acc);
% ACCmag25 = sqrt(curAcc_resampled(1,:).^2 + curAcc_resampled(2,:).^2 + curAcc_resampled(3,:).^2);
% ACCmag25_log(i) = mean(ACCmag25);
% Hacc(i) = entropy_proxy_hist(ACCmag25, nbits_entropy);
% if i == 1
%     dHacc(i) = 0;
% else
%     dHacc(i) = Hacc(i) - Hacc(i-1);
% end

% ACC control stream at fixed 125 Hz for entropy calculation
ACCmag25 = sqrt(curDataFilt2(1,:).^2 + curDataFilt2(2,:).^2 + curDataFilt2(3,:).^2);    %vector magnitude
accmag_val = mean(ACCmag25);   %temporal mean avg across window
hacc_val = entropy_proxy_hist(ACCmag25, nbits_entropy);
if i == 1
    dhacc_val = 0;
else
    dhacc_val = hacc_val - Hacc(i-1);
end

    %% Adaptive sampling controller logic (simple hysteresis on entropy jumps)

 % Decide on history to check
i0 = max(1, i - N_look_back + 1);
recent = [dHacc(i0:i-1), dhacc_val];  % Include current value
enough_hist = (i >= N_look_back);

% Default: no change
fs_next = fs_lo;  % Default to current mode's rate
next_state_mode = state_mode;
hi_timer_out = hi_timer_in;

% State machine controller
if FORCE_LOW
    state_mode = "LOW"; fs_next = fs_lo;
else
    switch state_mode
        case "LOW"
            % fs_adc = fs_lo;
            if enough_hist
                n_bad = sum(abs(recent) > Th_hi);         % count unstable windows in look-back
                if n_bad >= N_unstable
                    next_state_mode = "HIGH";
                    hi_timer_out = hi_hold;              % force HIGH for at least hi_hold windows
                    fs_next = fs_hi;
                end
            end
        case "HIGH"
            % fs_adc = fs_hi;
            % Decrement hold timer while in HIGH
            if hi_timer_out > 0
                hi_timer_out = hi_timer_out - 1;
            end

            % After hold expires, exit only if ALL last 5 are stable (< Th_low)
            if enough_hist && hi_timer_out == 0
                n_good = sum(abs(recent) < Th_low);   % count stable in last N_look_back
                if n_good >= N_stable
                    next_state_mode = "LOW";
                    fs_next = fs_lo;
                end
            end
        
    end
end

end