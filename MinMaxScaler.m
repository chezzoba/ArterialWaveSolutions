classdef MinMaxScaler
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        xmin;
        xmax;
    end

    methods
        function obj = MinMaxScaler(xmin, xmax)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if length(xmin) ~= length(xmax)
                error('Size of xmin has to be equal to that of xmax');
            end
            obj.xmin = xmin;
            obj.xmax = xmax;
        end

        function xprime = transform(obj, x)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(x) ~= length(obj.xmin)
                error('Size of x has to be equal to that of xmin & xmax');
            end
            xprime = (x - obj.xmin) ./ (obj.xmax - obj.xmin);
        end

        function x = inv_transform(obj, xp)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(xp) ~= length(obj.xmin)
                error('Size of xp has to be equal to that of xmin & xmax');
            end
            x = xp .* (obj.xmax - obj.xmin) + obj.xmin;
        end
    end
end