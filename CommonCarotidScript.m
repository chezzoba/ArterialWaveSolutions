params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
ccm = CommonCarotidModel;
zero = zeros(1, length(params));
[xmin, xmax, xact] = deal(zero, zero, zero);
for i = 1:length(params)
    xiv = ccm.(params(i));
    xmin(i) = xiv / 10;
    xmax(i) = xiv * 10;
    xact(i) = xiv;
end

scaler = MinMaxScaler(xmin, xmax);
load('Data/ccm_sol_model.mat');
% msol = ccm.model();
% spe = abs(sol - msol) .^ 2 / abs(sol) .^ 2;
% mspe = mean(spe, 2);
% err = sum(mspe);

erract = ccm.globalerr(scaler.transform(xact), scaler, params, sol);

SGD = GradientDescent;

fun = @(x) ccm.globalerr(x, scaler, params, sol);
x0 = 0.8 .* scaler.transform((xact));
dx = repmat([1e-5], length(xact));
options = optimset('PlotFcns',@optimplotfval);
xpred = scaler.inv_transform(fminsearch(fun, x0, options))
Perr = 100 .* abs(xpred - xact) ./ xact
%xpred = scaler.inv_transform(mean(SGD.optimize(fun, x0, dx), 1))

