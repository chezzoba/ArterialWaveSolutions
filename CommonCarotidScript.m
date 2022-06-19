addpath("Algorithms/", "Models/");

params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005,0.0001,0.002,0,0,0];
xmax = [0.3,0.02,0.005,740000000,5.22e+09,1e-7];
params = params(1:3);
xmin = xmin(1:3);
xmax = xmax(1:3);

scaler = MinMaxScaler(xmin, xmax);

problem = OptimisationProblem(CommonCarotidModel, params, scaler);


measurements = 0;

switch (measurements)
    case (0)
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.x0Tol = 1;
        problem.optimiser.epochs = 8;
        problem.model.optsol = [2, 3];
        
        [xpred, errP, erract, nguesses] = problem.fitsolution(0)
    case (1)
        problem.optimiser = ConOptimisation;
        problem.optimiser.x0Tol = 0.05;
        problem.optimiser.StepTolerance = 1e-4;
        problem.optimiser.TolFun = 2e-4;
        problem.model.optsol = [2, 3];

        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end