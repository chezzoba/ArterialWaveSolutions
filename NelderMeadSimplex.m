classdef NelderMeadSimplex

    properties
        TolX = 1e-10;
        TolFun = 1e-21;
        epochs = 5;
        plt = true;
        x0Tol = 4;
        lenX = 6;
    end

    methods

        function xpred = fit(obj, fun)
            %Fit - Optimise function fun, with an initial guess x0
            if (obj.plt)
                options = optimset('PlotFcns',@optimplotfval,...
                'TolX', obj.TolX, 'TolFun', obj.TolFun);
            else
                options = optimset('TolX', obj.TolX, 'TolFun', obj.TolFun);
            end
            
            x0 = rand(1, obj.lenX);
            nguesses = 1;
            while fun(x0) > obj.x0Tol
                x0 = rand(1, obj.lenX);
                nguesses = nguesses + 1;
            end
            nguesses

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