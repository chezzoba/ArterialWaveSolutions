addpath("Algorithms/");

Z = P11out ./ Q11out;
om = omegas;

zfun = @(x) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
    ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]


opt = NelderMeadSimplex;
scaler = MinMaxScaler([0, 0, 0], [1e11,1e11,1e-7]);

xact = [RW1(11), RW2(11), CWK(11)];

fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp)) - Z).^2);
opt.lenX = 3;
opt.x0Tol = 5e16;
xpred = opt.fit(fun, 0);
Perr = 100.*opt.relerr(scaler.inv_transform(xpred), xact)