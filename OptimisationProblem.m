classdef OptimisationProblem
    %OPTIMISATIONPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model;
        params;
        scaler;
        xact;
        optimiser = ConOptimisation;
    end
    
    methods
        function obj = OptimisationProblem(model, params, scaler)
            %OPTIMISATIONPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            [obj.model, obj.params, obj.scaler] = ...
                deal(model, params, scaler);
            obj.xact = zeros(1, length(params));
            for i = 1:length(params)
                obj.xact(i) = obj.model.(params(i));
            end
        end
        
        function [xpred, errP, erract, nguesses] = fitmeasurements(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xactt = obj.scaler.transform(obj.xact);
            erract = obj.model.measurementerr(xactt, obj.scaler, obj.params);
            fun = @(x) obj.model.measurementerr(x, obj.scaler, obj.params);
            [xpred, nguesses] = obj.optimiser.fit(fun, 0);
            xpredt = obj.scaler.inv_transform(xpred);
            errP = 100 .* obj.optimiser.relerr(obj.xact, xpredt);
        end

        function [xpred, errP, erract, nguesses] = fitsolution(obj, sol)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xactt = obj.scaler.transform(obj.xact);
            erract = obj.model.globalerr(xactt, obj.scaler, obj.params, sol);
            fun = @(x) obj.model.globalerr(x, obj.scaler, obj.params, sol);
            [xpred, nguesses] = obj.optimiser.fit(fun, 0);
            xpredt = obj.scaler.inv_transform(xpred);
            errP = 100 .* obj.optimiser.relerr(obj.xact, xpredt);
        end
    end
end

