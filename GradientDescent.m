classdef GradientDescent
    %GradientDescent Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        max_iter = 1e4; % Maximum number of iterations
        tol = 1e-3; % Tolerance
        alpha = 1e-2; % Step Size
    end
    
    methods
        function x1 = optimize(obj, fun, x0, dx)
            wbar = waitbar(0,'Loading Optimisation...','Name','Optimisation Running',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            setappdata(wbar,'canceling',0);

            for iteration = 1:obj.max_iter
                if getappdata(wbar,'canceling')
                    break;
                end
                [g, funh] = GradientDescent.grad(fun, x0, dx);
                a = obj.alpha ./ norm(g);
                x1 = x0 - a .* g;
                if all(abs(x1 - x0) < obj.tol)
                    break;
                else
                    x0 = x1;
                end
                waitbar(iteration / obj.max_iter, wbar, sprintf('Error: %12.9f', funh))
            end
            delete(wbar);
        end
    end

    methods (Static)
        function [g, funh] = grad(fun, x, dx)
            g = zeros(1, length(x));
            for i = 1:length(x)
                [xh, xl] = deal(x, x);
                [xh(i), xl(i)] = deal(x(i) + dx(i) / 2, x(i) - dx(i) / 2);
                funh = fun(xh);
                g(i) = (funh - fun(xl)) / dx(i);
            end
        end
    end
end

