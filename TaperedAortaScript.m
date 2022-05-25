addpath("Algorithms/");
addpath("Models/");

params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005,0,0.002,0,0,0];
xmax = [0.3,0.0152,0.005,740000000,5.22e+09,1e-7];
scaler = MinMaxScaler(xmin, xmax);

problem = OptimisationProblem(TaperedAortaModel, params, scaler);

[xpred, errP, erract, nguesses] = problem.fitmeasurements(0)

% problem.optimiser = NelderMeadSimplex;
% problem.optimiser.x0Tol = 8;

% [xpred, errP, erract, nguesses] = problem.fitsolution(problem.model.model(), 0)