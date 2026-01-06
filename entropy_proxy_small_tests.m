%% Simple sanity checks for entropy proxy helpers
% Run: >> entropy_proxy_small_tests

rng(0);           % reproducible noise
nbits = 3;

constSig   = zeros(1, 64); % all zeros
stepSig    = [zeros(1, 32) ones(1, 32)]; % half zeros, half ones
periodicSig = repmat([0 1 0 1 2 1], 1, 11); % short repeating pattern
noisySig   = randn(1, 64); % white Gaussian noise

proxies = {@entropy_proxy_hist, @entropy_proxy_arith_o1, @entropy_proxy_lz};
names   = {'hist', 'arith_o1', 'lz'};

margin = 1e-5;    % ordering tolerance
invarTol = 1e-2;  % scale/offset invariance tolerance

for k = 1:numel(proxies)
    f = proxies{k};
    name = names{k};

    H_const   = f(constSig, nbits);
    H_step    = f(stepSig, nbits);
    H_step_sh = f(stepSig + 5, nbits);  % shifted version should match
    H_period  = f(periodicSig, nbits);
    H_noisy   = f(noisySig, nbits);

    fprintf('%s: const=%.4f, step=%.4f, step+offset=%.4f, periodic=%.4f, noisy=%.4f\n', ...
        name, H_const, H_step, H_step_sh, H_period, H_noisy);

    assert(H_noisy > H_step - margin,     '%s: noisy should be higher than step', name);
    assert(H_noisy > H_period - margin,   '%s: noisy should be higher than periodic', name);
    assert(H_noisy > H_const - margin,    '%s: noisy should be higher than constant', name);
    assert(H_step  >= H_const - margin,   '%s: step should be higher than constant', name);
    assert(abs(H_step - H_step_sh) < invarTol, '%s: scale/offset invariance failed', name); % offset should not change entropy
end

fprintf('All sanity checks passed.\n');
