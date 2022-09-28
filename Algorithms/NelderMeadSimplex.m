classdef NelderMeadSimplex
% NelderMeadSimplex Unconstrained Optimisation Algorithm using
% Nelder Mead Simplex Search algorithm
    properties
        TolX = 1e-50; % Minimum change in x allowed before breaking
        TolFun = 1e-50; % Minimum change in function allowed before breaking
        epochs = 3; % Number of optimisations performed in sequence
        plt = true; % Plot function value vs iteration
        x0Tol = 5; % Constraint of the suitability of the guess
        lenX = 6; % Length of input vector
        MaxFunEvals = 1000; % Maximum number of function evaluations allowed
        MaxIter = 1000; % Maximum number of iterations allowed
    end

    methods

        function [xpred, nguesses] = fit(obj, fun, x0)
            %Fit - Optimise function fun, with an initial guess x0
            % If x0 invalid the algorithm will try to guess x0 based on
            % the value set in x0Tol
            if (obj.plt)
                options = optimset('PlotFcns',@optimplotfval,...
                'TolX', obj.TolX, 'TolFun', obj.TolFun,...
                'MaxFunEvals', obj.MaxFunEvals, 'MaxIter', obj.MaxIter);
            else
                options = optimset('TolX', obj.TolX, 'TolFun', obj.TolFun,...
                    'MaxFunEvals', obj.MaxFunEvals, 'MaxIter', obj.MaxIter);
            end

            nguesses = 0;
            if (length(x0) < obj.lenX)
                x0 = rand(1, obj.lenX);
                while fun(x0) > obj.x0Tol
                    x0 = rand(1, obj.lenX);
                    nguesses = nguesses + 1;
                end
            end


            xpred = x0;
            for epoch = 1:obj.epochs
                xpred = fminsearch(fun, xpred, options);
            end
        end
    end

    methods (Static)

        function rel_error = relerr(Pact, Ppred)
            rel_error = abs((Pact - Ppred) ./ Pact);
        end

    end
end