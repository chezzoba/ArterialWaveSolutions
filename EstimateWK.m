addpath("Algorithms/");

RW1s = [0,0,0,0,0,0,0,0,0,0,5.1918,19.1515,9.882,11.7617,17.4352,34.1378,34.1378,74.0167,5.9149,5.9149]*(10^7);
RW2s = [0,0,0,0,0,0,0,0,0,0,10.6080,52.2129,13.0183,7.5726,5.5097,5.3949,5.3949,46.2252,10.1737,10.1737]*(10^8);
CWKs = [0,0,0,0,0,0,0,0,0,0,8.6974,1.767,7.0871,12.1836,16.7453,17.1017,17.1017,1.9959,9.0686,9.0686]*(10^-10);

lsub_outlet_13_flow = load('PreviousCode/automeris/aorta/lsub_outlet_13_flow.csv');
lsub_outlet_13_Pressure = load('PreviousCode/automeris/aorta/lsub_outlet_13_Pressure.csv');

NP = length(lsub_outlet_13_Pressure);
Pout = fft(lsub_outlet_13_Pressure(:, 2));
omP = 2*pi/(lsub_outlet_13_Pressure(end, 1)) * (0:NP);

NQ = length(lsub_outlet_13_flow);
Qout = fft(lsub_outlet_13_flow(:, 2));
omQ = 2*pi/(lsub_outlet_13_flow(end, 1)) * (0:NQ);

om = omP(3:end-1);

Pout = interp1(omP(2:end-1), Pout(2:length(omP)-1), om);
Qout = interp1(omQ(2:end-1), Qout(2:length(omQ)-1), om);

Z = Pout ./ Qout;

zfun = @(x) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
    ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]


opt = NelderMeadSimplex;
scaler = MinMaxScaler([0, 0, 0], [1e11,1e11,1e-7]);

xact = [RW1s(13), RW2s(13), CWKs(13)];

fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp)) - Z).^2);
opt.lenX = 3;
opt.x0Tol = 5e16;
xpred = opt.fit(fun, xact);
Perr = 100.*opt.relerr(scaler.inv_transform(xpred), xact)