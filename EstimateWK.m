addpath("Algorithms/");

Z = Z5;
om = omegas;

zfun = @(x, om) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
    ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]


opt = NelderMeadSimplex;
scaler = MinMaxScaler([0, 0, 0], [740000000,5.22e+09,1e-7]);

fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z).^2);
opt.lenX = 3;
opt.x0Tol = 1e9;
[xpred, nguesses] = opt.fit(fun, [0.5, 0.5, 0.5]);