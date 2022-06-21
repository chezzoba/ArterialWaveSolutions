classdef SingleVessel
    %SINGLEVESSEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R;
        L;
        a;
        beta;
        rho;
    end
    
    methods
        function obj = SingleVessel(R, L, a, beta, rho)
            %SINGLEVESSEL Construct an instance of this class
            %   Detailed explanation goes here
            obj.R, obj.L, obj.a, obj.beta, obj.rho = ...
                deal(R, L, a, beta, rho);
        end
        
        function sv = s(obj, x)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            sv = (obj.R - tan(obj.a)*x)/sin(obj.a);
        end

        function [B_A, Yeff] = backpropagate(obj, omega)
            s4out = obj.s(obj.L);
            [J_13s4out,Y_13s4out,J_43s4out,Y_43s4out,fs4out] = besselfunctions(obj.a,s4out,omega,obj.rho,obj.beta);    
            Y_s4out = (2*pi*(1-cos(obj.a)))*(fs4out/obj.rho)^0.5*s4out^2.5;
            B_A = -(1i*Y_s4out*J_43s4out+Yeff_5*J_13s4out)/(Yeff_5*Y_13s4out+1i*Y_s4out*Y_43s4out);
            
            s4in = obj.s(0);
            [J_13s4in,Y_13s4in,J_43s4in,Y_43s4in,fs4in] = besselfunctions(obj.a,s4in,omega,obj.rho,obj.beta);      
            Y_s4in = (2*pi*(1-cos(obj.a)))*(fs4in/obj.rho)^0.5*s4in^2.5;
            Yeff = -1i*Y_s4in*(J_43s4in+B_A*Y_43s4in)/(J_13s4in+B_A*Y_13s4in);
        end

        function [] = forwardpropagate(obj, omega)
            
        end
    end
end

