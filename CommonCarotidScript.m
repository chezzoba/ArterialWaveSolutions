%% Common Carotid Artery Parameter Estimation
addpath("Algorithms/", "Models/");

% Setting Parameters to be estimated and physical ranges thereof
params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005,0.0001,0.002,0,0,0];
xmax = [0.3,0.02,0.005,740000000,5.22e+09,1e-7];
scaler = MinMaxScaler(xmin, xmax);

% Defining estimation problem
problem = OptimisationProblem(CommonCarotidModel, params, scaler);

% Option to choose what values to fit the model to
% 0 - fit to 1D model output using known inputs
% 1 - fit to 3D numerical simulation data
measurements = 0;

switch (measurements)
    case (0)
        % Choosing optimisation algorithm and tuning parameters
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.x0Tol = 0.001;
        problem.optimiser.epochs = 10;
        problem.model.optsol = [2];
        
        % Fit model to 1D model outputs using known inputs
        [xpred, errP, erract, nguesses] = problem.fitsolution(0)
    case (1)
        % Choosing optimisation algorithm and tuning parameters
        problem.optimiser = ConOptimisation;
        problem.optimiser.x0Tol = 0.05;
        problem.optimiser.StepTolerance = 1e-4;
        problem.optimiser.TolFun = 2e-4;
        problem.model.optsol = [2, 3];
        
        % Fit model to 3D numerical simulation data
        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end