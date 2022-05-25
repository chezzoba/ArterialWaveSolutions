classdef NelderMeadSimplex

    properties
        TolX = 1e-22;
        TolFun = 2e-16;
        epochs = 1;
        plt = true;
        x0Tol = 0.00001;
        lenX = 6;
    end

    methods

        function [xpred, nguesses] = fit(obj, fun, x0)
            %Fit - Optimise function fun, with an initial guess x0
            if (obj.plt)
                options = optimset('PlotFcns',@optimplotfval,...
                'TolX', obj.TolX, 'TolFun', obj.TolFun);
            else
                options = optimset('TolX', obj.TolX, 'TolFun', obj.TolFun);
            end

            nguesses = 1;
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