addpath("Algorithms/", "Models/", "VesselModels/");

params = ["RW17" "RW27" "Cwk7" "RW18" "RW28" "Cwk8" "RW19" "RW29" "Cwk9" "RW110" "RW210" "Cwk10"];
xmin = repmat([0, 0, 0], 1, 4);
xmax = repmat([1e12,1e12,1e-6], 1, 4);
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
        problem.optimiser.x0Tol = 1e4;

        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end

xpredt = scaler.inv_transform(xpred);

