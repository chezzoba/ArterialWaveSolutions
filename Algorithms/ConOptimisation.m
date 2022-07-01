classdef ConOptimisation
    %CONOPTIMISATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        plt = true;
        lenX = 6;
        x0Tol = 1;
        TolX = 1e-30;
        TolFun = 1e-5;
        StepTolerance = 1e-50;
    end
    
    methods
        function [xpred, nguesses] = fit(obj, fun, x0)
            if (obj.plt)
                options = optimoptions(@fmincon, 'PlotFcns', @optimplotfval,...
    'TolX', obj.TolX, 'TolFun', obj.TolFun, 'StepTolerance', obj.StepTolerance,...
    'OptimalityTolerance', 1e-20, 'MaxFunctionEvaluations', 3e3*length(obj.lenX));
            else
                options = optimoptions(@fmincon, 'TolX', obj.TolX,...
                    'TolFun', obj.TolFun, 'StepTolerance', obj.StepTolerance,...
                    'OptimalityTolerance', 1e-20, 'MaxFunctionEvaluations', 3e3*length(obj.lenX));
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

