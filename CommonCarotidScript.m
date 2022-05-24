params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005, 3e-3, 0.0013, 0, 0, 0];
xmax = [0.15, 15.2e-3, 0.005, 74e7, 52.2e8, 17.1e-10];

scaler = MinMaxScaler(xmin, xmax);
load('Data/ccm_sol_model.mat');

optimise = OptimisationProblem(CommonCarotidModel, params, scaler);

[xpred, errP, erract, nguesses] = optimise.fitmeasurements()

% optimise.optimiser.x0Tol = 3;
% optimise.optimiser.TolFun = 1e-25;
% optimise.optimiser.StepTolerance = 0;
% 
% [xpred, errP, erract, nguesses] = optimise.fitsolution(sol)