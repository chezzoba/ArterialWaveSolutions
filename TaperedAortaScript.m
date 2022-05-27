addpath("Algorithms/", "Models/");

params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005,0,0,0,0,0];
xmax = [0.3,0.02,0.005,740000000,5.22e+09,1e-7];
scaler = MinMaxScaler(xmin, xmax);

problem = OptimisationProblem(TaperedAortaModel, params, scaler);


measurements = 1;

switch (measurements)
    case (0)
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.x0Tol = 5;
        problem.optimiser.epochs = 8;
        
        [xpred, errP, erract, nguesses] = problem.fitsolution(0)
    case (1)
        problem.optimiser = ConOptimisation;
        problem.optimiser.x0Tol = 1;

        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end

