addpath("Algorithms/");
addpath("Models/");

params = ["L" "R" "be" "RW1" "RW2" "Cwk"];
xmin = [0.005,0,0,0,0,0];
xmax = [0.3,0.02,0.005,740000000,5.22e+09,1e-7];
scaler = MinMaxScaler(xmin, xmax);
model = TaperedCarotidModel;
problem = OptimisationProblem(TaperedCarotidModel, params, scaler);


measurements = 1;

switch (measurements)
    case (0)
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.epochs = 3;
        problem.optimiser.TolFun = 2e-29;
        problem.optimiser.x0Tol = 1;
        
        [xpred, errP, erract, nguesses] = problem.fitsolution(0)
    case (1)
        problem.optimiser = ConOptimisation;
        problem.optimiser.x0Tol = 1;
        
        [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)
end

