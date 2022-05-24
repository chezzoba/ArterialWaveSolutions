params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
ccm = CommonCarotidModel;

xmin = [0.005, 3e-3, 0.0013, 0, 0, 0];
xmax = [0.15, 15.2e-3, 0.005, 74e7, 52.2e8, 17.1e-10];

scaler = MinMaxScaler(xmin, xmax);
load('Data/ccm_sol_model.mat');

erract = ccm.globalerr(scaler.transform(xact), scaler, params, sol);

fun = @(x) ccm.globalerr(x, scaler, params, sol);


nms = NelderMeadSimplex;

[xpred, nguesses] = nms.fit(fun)

errP = 100 .* nms.relerr(xact, scaler.inv_transform(xpred))

