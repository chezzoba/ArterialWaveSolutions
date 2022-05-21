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

scaler = MinMaxScaler(xmin, xmax);
load('Data/ccm_sol_model.mat');
% msol = ccm.model();
% spe = abs(sol - msol) .^ 2 / abs(sol) .^ 2;
% mspe = mean(spe, 2);
% err = sum(mspe);

erract = ccm.globalerr(scaler.transform(xact), scaler, params, sol);


fun = @(x) ccm.globalerr(x, scaler, params, sol);
x00 = scaler.transform(2 .* (xact));
options = optimset('PlotFcns',@optimplotfval, 'TolX', 1e-10, 'TolFun', 1e-21);

xpred = fminsearch(fun, x00, options);

total_iters = 10;

for iter = 1:total_iters
    x0 = xpred;
    xpred = fminsearch(fun, x0, options);
end

err0 = 100 .* abs(scaler.inv_transform(x00) - xact) ./ xact
errP = 100 .* abs(scaler.inv_transform(xpred) - xact) ./ xact

