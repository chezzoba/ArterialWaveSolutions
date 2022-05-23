params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
ccm = CommonCarotidModel;
% zero = zeros(1, length(params));
% [xmin, xmax, xact] = deal(zero, zero, zero);
% for i = 1:length(params)
%     xiv = ccm.(params(i));
%     xmin(i) = xiv / 5;
%     xmax(i) = xiv * 5;
%     xact(i) = xiv;
% end
xmin = [0.005, 3e-3, 0.0013, 0, 0, 0];
xmax = [0.15, 15.2e-3, 0.005, 74e7, 52.2e8, 17.1e-10];

scaler = MinMaxScaler(xmin, xmax);
load('Data/ccm_sol_model.mat');
% msol = ccm.model();
% spe = abs(sol - msol) .^ 2 / abs(sol) .^ 2;
% mspe = mean(spe, 2);
% err = sum(mspe);

erract = ccm.globalerr(scaler.transform(xact), scaler, params, sol);

fun = @(x) ccm.globalerr(x, scaler, params, sol);


nms = NelderMeadSimplex;

xpred = nms.fit(fun)

errP = 100 .* nms.relerr(xact, scaler.inv_transform(xpred))

