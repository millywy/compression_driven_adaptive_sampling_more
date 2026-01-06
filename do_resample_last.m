function sig_res = do_resample_last(sig, fs_in, fs_out)
%DO_RESAMPLE Resample multi-channel signal (channels x samples) with proper anti-aliasing.
% Uses DECIMATE for integer downsampling ratios, otherwise RESAMPLE with rational p/q.
%
% sig:   [nCh x nSamp] (time along columns)
% fs_in: input sample rate (Hz)
% fs_out: output sample rate (Hz)

    if abs(fs_out - fs_in) < 1e-12
        sig_res = sig;
        return;
    end

    [nCh, nSamp] = size(sig);

    % ---------- Case 1: integer-factor downsample ----------
    r = fs_in / fs_out;
    if fs_out < fs_in && abs(r - round(r)) < 1e-9
        decim = round(r);
        % decimate includes an anti-alias lowpass; FIR is safer for bio signals
        outLens = ceil(nSamp / decim);
        sig_res = zeros(nCh, outLens);
        for c = 1:nCh
            y = decimate(sig(c,:), decim, 'fir');
            sig_res(c, 1:numel(y)) = y;
        end
        % trim in case decimate returned slightly longer
        sig_res = sig_res(:, 1:max(1, min(size(sig_res,2), outLens)));
        return;
    end

    % ---------- Case 2: integer-factor upsample ----------
    u = fs_out / fs_in;
    if fs_out > fs_in && abs(u - round(u)) < 1e-9
        up = round(u);
        % resample has a good interpolation/anti-imaging filter
        sig_res = resample(sig.', up, 1).';
        return;
    end

    % ---------- Case 3: rational resampling ----------
    ratio = fs_out / fs_in;
    [p,q] = rat(ratio, 1e-12);
    sig_res = resample(sig.', p, q).';
end
