addpath("Algorithms/");
addpath("Models/");
load("Problems/CommonCarotidArtery.mat");

% [xpred, errP, erract, nguesses] = problem.fitmeasurements(0)


problem.optimiser = NelderMeadSimplex;
problem.optimiser.x0Tol = 1;

[xpred, errP, erract, nguesses] = problem.fitsolution(problem.model.model(),0)