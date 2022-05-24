params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
ccm = CommonCarotidModel;
zero = zeros(1, length(params));
[xmin, xmax, xact] = deal(zero, zero, zero);
for i = 1:length(params)
    xiv = ccm.(params(i));
    xmin(i) = xiv / 5;
    xmax(i) = xiv * 5;
    xact(i) = xiv;
end


xmin = [0.005, 3e-3, 0.0013, 0, 0, 0];
xmax = [0.15, 15.2e-3, 0.005, 74e7, 52.2e8, 17.1e-10];

scaler = MinMaxScaler(xmin, xmax);
%load('Data/ccm_sol_model.mat');
sol = ccm.measurement();

%erract = ccm.globalerr(scaler.transform(xact), scaler, params, sol)
erract = ccm.measurementerr(scaler.transform(xact), scaler, params)

%fun = @(x) ccm.globalerr(x, scaler, params, sol);
fun = @(x) ccm.measurementerr(x, scaler, params);


nms = NelderMeadSimplex;

nguesses = 1;

x0 = rand(1, length(xact));
while fun(x0) > 5e-2
    x0 = rand(1, length(xact));
    nguesses = nguesses + 1;
end

options = optimoptions(@fmincon, 'PlotFcns', @optimplotfval,...
    'TolX', 1e-30, 'TolFun', 1e-5, 'StepTolerance', 1e-20);

problem = struct('x0', x0, 'objective', fun, 'lb', zeros(1, length(x0)),...
    'ub', ones(1, length(x0)), 'solver', 'fmincon', 'options', options);

xpred = fmincon(problem)
xpredt = scaler.inv_transform(xpred);

% [xpred, nguesses] = nms.fit(fun, 0)

errP = 100 .* nms.relerr(xact, xpredt)

