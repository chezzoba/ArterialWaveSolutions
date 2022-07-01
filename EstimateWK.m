addpath("Algorithms/");

RW1s = [0,0,0,0,0,0,0,0,0,0,5.1918,19.1515,9.882,11.7617,17.4352,34.1378,34.1378,74.0167,5.9149,5.9149]*(10^7);
RW2s = [0,0,0,0,0,0,0,0,0,0,10.6080,52.2129,13.0183,7.5726,5.5097,5.3949,5.3949,46.2252,10.1737,10.1737]*(10^8);
CWKs = [0,0,0,0,0,0,0,0,0,0,8.6974,1.767,7.0871,12.1836,16.7453,17.1017,17.1017,1.9959,9.0686,9.0686]*(10^-10);

flows = load('PreviousCode/automeris/aorta/bc_outlet_11_flow.csv');
pressures = load('PreviousCode/automeris/aorta/bc_outlet_11_Pressure.csv');

Z = P7out ./ Q7out;
om = omegas;

zfun = @(x) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
    ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]


opt = ConOptimisation;
opt.plt = false;
scaler = MinMaxScaler([1e7, 1e8, 0], [1e9,1e10,1e-8]);

fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp)) - Z).^2);
opt.lenX = 3;
opt.x0Tol = 1e30;
xpred = scaler.inv_transform(opt.fit(fun, 0))

% xact = [RW1s(11), RW2s(11), CWKs(11)]

% Perr = 100.*ConOptimisation.relerr(xpred, xact)