addpath("Algorithms/", "Models/", "VesselModels/");

params = ["RW17" "RW27" "Cwk7"];
xmin = [0,0,0];
xmax = [1e11,1e11,1e-6];
scaler = MinMaxScaler(xmin, xmax);

problem = OptimisationProblem(AorticArchModel, params, scaler);


measurements = 1;

switch (measurements)
    case (0)
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.x0Tol = 5;
        problem.optimiser.epochs = 8;
        
        [xpred, errP, erract, nguesses] = problem.fitsolution(0)
    case (1)
        problem.optimiser = ConOptimisation;
        problem.optimiser.x0Tol = 5;

        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end

