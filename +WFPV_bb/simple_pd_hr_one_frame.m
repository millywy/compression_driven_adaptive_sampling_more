function [BPM_val, state] = simple_pd_hr_one_frame(curDataRaw, fs0, fs_adc, fs_proc, FFTres, WFlength, searchHz, state, i, BPM_est, idnb, bpHz)
% Simple HR estimator: bandpass -> resample -> peak detection in PPG.

% Unused inputs kept for drop-in compatibility.
% FFTres, WFlength, idnb are not used in this simple estimator.

% if ~isfield(state, 'prevBPM')
%     state.prevBPM = NaN;
% end

% Bandpass filter PPG at original rate
[b, a] = butter(4, bpHz / (fs0 / 2), 'bandpass');
PPG_raw = curDataRaw(1:2, :);
PPG_filt = zeros(size(PPG_raw));
for c = 1:size(PPG_raw, 1)
    PPG_filt(c, :) = filter(b, a, PPG_raw(c, :));
end

% Resample to target ADC rate
PPG_adc = do_resample_last(PPG_filt, fs0, fs_adc);

% Optional internal processing rate
if abs(fs_adc - fs_proc) < eps
    PPG = PPG_adc;
    fs = fs_adc;
else
    PPG = do_resample_last(PPG_adc, fs_adc, fs_proc);
    fs = fs_proc;
end

% Normalize and blend the two PPG channels
PPG1 = PPG(1, :);
PPG2 = PPG(2, :);
PPG1 = (PPG1 - mean(PPG1)) / (std(PPG1) + eps);
PPG2 = (PPG2 - mean(PPG2)) / (std(PPG2) + eps);
PPG_ave = 0.5 * (PPG1 + PPG2);

% Light smoothing for peak detection stability
win = max(1, round(0.15 * fs));
if win > 1
    PPG_ave = filter(ones(1, win) / win, 1, PPG_ave);
    PPG_ave = (PPG_ave - mean(PPG_ave)) / (std(PPG_ave) + eps);
end

% Peak detection
minDist = max(1, floor(fs / max(searchHz(2), eps)));
minProm = 0.3; % normalized units
[~, locs] = findpeaks(PPG_ave, 'MinPeakDistance', minDist, 'MinPeakProminence', minProm);

% % If needed, try inverted signal and keep the richer set
% if numel(locs) < 2
%     [~, locs_inv] = findpeaks(-PPG_ave, 'MinPeakDistance', minDist, 'MinPeakProminence', minProm);
%     if numel(locs_inv) > numel(locs)
%         locs = locs_inv;
%     end
% end

% Calculate HR from inter-beat intervals
BPM_val = NaN;
if numel(locs) >= 2
    ibi = diff(locs) / fs; % in seconds
    minIbi = 1 / searchHz(2);
    maxIbi = 1 / searchHz(1);
    ibi = ibi(ibi >= minIbi & ibi <= maxIbi);
    if ~isempty(ibi)
        BPM_val = 60 / median(ibi);
    else
        duration = (locs(end) - locs(1)) / fs;
        if duration > 0
            BPM_val = 60 * (numel(locs) - 1) / duration;
        end
    end
end

% Keep within the search band
if ~isnan(BPM_val)
    minBPM = 60 * searchHz(1);
    maxBPM = 60 * searchHz(2);
    if BPM_val < minBPM || BPM_val > maxBPM
        BPM_val = NaN;
    end
end

% % Fallback to previous estimate if needed
% if isnan(BPM_val)
%     if ~isnan(state.prevBPM)
%         BPM_val = state.prevBPM;
%     elseif i > 1 && numel(BPM_est) >= i - 1 && ~isnan(BPM_est(i - 1)) && BPM_est(i - 1) > 0
%         BPM_val = BPM_est(i - 1);
%     end
% end

% if ~isnan(BPM_val)
%     state.prevBPM = BPM_val;
% end

end
