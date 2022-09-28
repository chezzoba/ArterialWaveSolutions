classdef OptimisationProblem
    %OPTIMISATIONPROBLEM Defines an optimisation problem
    % Where the magic happens
    
    properties
        model; % Model (from Models/) used in the estimation technique.
        params; % Array of names of parameters to be optimised
        scaler; % MinMaxScaler object corresponding to params
        xact; % Known parameter values
        optimiser = ConOptimisation; % Optimiser used
    end
    
    methods
        function obj = OptimisationProblem(model, params, scaler)
            %OPTIMISATIONPROBLEM Initialise the problem
            [obj.model, obj.params, obj.scaler] = deal(model, params, scaler);
            obj.xact = zeros(1, length(params));
            for i = 1:length(params)
                obj.xact(i) = obj.model.(params(i));
            end
        end
        
        function [xpred, errP, erract, nguesses] = fitmeasurements(obj, x0)
            % Fits model to 3D numerical simulation measurements
            obj.optimiser.lenX = length(obj.params);
            fun = @(x) obj.model.measurementerr(x, obj.scaler, obj.params);
            erract = fun(obj.scaler.transform(obj.xact));
            [xpred, nguesses] = obj.optimiser.fit(fun, x0);
            xpredt = obj.scaler.inv_transform(xpred);
            errP = 100 .* obj.optimiser.relerr(obj.xact, xpredt);
        end

        function [xpred, errP, erract, nguesses] = fitsolution(obj, x0)
            % Fits model to model solutions that use the known parameters
            sol = obj.model.model();
            obj.optimiser.lenX = length(obj.params);
            fun = @(x) obj.model.globalerr(x, obj.scaler, obj.params, sol);
            erract = fun(obj.scaler.transform(obj.xact));
            [xpred, nguesses] = obj.optimiser.fit(fun, x0);
            xpredt = obj.scaler.inv_transform(xpred);
            errP = 100 .* obj.optimiser.relerr(obj.xact, xpredt);
        end
    end
end

