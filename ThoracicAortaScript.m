addpath("Algorithms/", "Models/");

params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005,0,0,0,0,0];
xmax = [0.3,0.02,0.0025,740000000,5.22e+09,1e-7];
scaler = MinMaxScaler(xmin, xmax);

problem = OptimisationProblem(ThoracicAortaModel, params, scaler);


measurements = 1;

switch (measurements)
    case (0)
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.x0Tol = 5;
        problem.optimiser.epochs = 10;
        
        [xpred, errP, erract, nguesses] = problem.fitsolution(0)
    case (1)
        problem.optimiser = ConOptimisation;
        problem.optimiser.x0Tol = 0.5;
        problem.optimiser.StepTolerance = 4e-4;
        problem.optimiser.TolFun = 2e-3;

        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end

