classdef NelderMeadSimplex

    properties
        TolX = 1e-50;
        TolFun = 1e-50;
        epochs = 3;
        plt = true;
        x0Tol = 5;
        lenX = 6;
        MaxFunEvals = 1000;
        MaxIter = 1000;
    end

    methods

        function [xpred, nguesses] = fit(obj, fun, x0)
            %Fit - Optimise function fun, with an initial guess x0
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