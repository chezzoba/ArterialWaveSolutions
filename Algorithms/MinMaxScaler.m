classdef MinMaxScaler
    % Scaler that normalises any vector based in its minimum and maximum
    % value

    properties
        xmin;
        xmax;
    end

    methods
        function obj = MinMaxScaler(xmin, xmax)
            %Initialise the scaler with min and max vector
            if length(xmin) ~= length(xmax)
                error('Size of xmin has to be equal to that of xmax');
            end
            obj.xmin = xmin;
            obj.xmax = xmax;
        end

        function xprime = transform(obj, x)
            % transform x to 0 and 1
            if length(x) ~= length(obj.xmin)
                error('Size of x has to be equal to that of xmin & xmax');
            end
            xprime = (x - obj.xmin) ./ (obj.xmax - obj.xmin);
        end

        function x = inv_transform(obj, xp)
            %Inverse transform xp from 0 and 1
            if length(xp) ~= length(obj.xmin)
                error('Size of xp has to be equal to that of xmin & xmax');
            end
            x = xp .* (obj.xmax - obj.xmin) + obj.xmin;
        end
    end
end