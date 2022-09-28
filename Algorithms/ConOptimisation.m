classdef ConOptimisation
    %CONOPTIMISATION Constrained Optimisation Algorithm using
    % the interior point method
    properties
        plt = true; % Plot function value vs iteration
        lenX = 6; % Length of input vector
        x0Tol = 1; % Constraint of the suitability of the guess
        TolX = 1e-30; % Minimum change in x allowed before breaking
        TolFun = 1e-5; % Minimum change in function allowed before breaking
        StepTolerance = 1e-50; % Minimum step size allowed before breaking
    end
    
    methods
        function [xpred, nguesses] = fit(obj, fun, x0)
            % Optimises a function, fun with initial guess x0.
            % If x0 invalid the algorithm will try to guess x0 based on
            % the value set in x0Tol
            if (obj.plt)
                options = optimoptions(@fmincon, 'PlotFcns', @optimplotfval,...
    'TolX', obj.TolX, 'TolFun', obj.TolFun, 'StepTolerance', obj.StepTolerance,...
    'OptimalityTolerance', 1e-20, 'MaxFunctionEvaluations', 6e3*length(obj.lenX));
            else
                options = optimoptions(@fmincon, 'TolX', obj.TolX,...
                    'TolFun', obj.TolFun, 'StepTolerance', obj.StepTolerance,...
                    'OptimalityTolerance', 1e-20, 'MaxFunctionEvaluations', 6e3*length(obj.lenX));
            end

            nguesses = 0;
            if (length(x0) < obj.lenX)
                x0 = rand(1, obj.lenX);
                nguesses = 1;
                while fun(x0) > obj.x0Tol
                    x0 = rand(1, obj.lenX);
                    nguesses = nguesses + 1;
                end
            end

        
        problem = struct('x0', x0, 'objective', fun, 'lb', zeros(1, length(x0)),...
            'ub', ones(1, length(x0)), 'solver', 'fmincon', 'options', options);
        xpred = fmincon(problem);
        end
    end

    methods (Static)

        function rel_error = relerr(Pact, Ppred)
            rel_error = abs((Pact - Ppred) ./ Pact);
        end

    end
end

