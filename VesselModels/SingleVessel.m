classdef SingleVessel
    %SINGLEVESSEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = 3; % 2, 3 or 5
        R;
        L;
        a;
        beta;
        rho;
        RW1;
        RW2;
        Cwk;
    end
    
    methods
        function obj = SingleVessel(R, L, a, beta, rho, WKP)
            %SINGLEVESSEL Construct an instance of this class
            %   Detailed explanation goes here
            obj.type = type;
            [obj.R, obj.L, obj.a, obj.beta, obj.rho] = ...
                deal(R, L, a, beta, rho);
            [obj.RW1, obj.RW2, obj.Cwk] = deal(WKP(1), WKP(2), WKP(3));

        end
        
        function sv = s(obj, x)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            sv = (obj.R - tan(obj.a)*x)/sin(obj.a);
        end

        function fv = f(obj)
            fv = obj.beta * tan(obj.a)^2 * sin(obj.a) / (1 - cos(obj.a));
        end

        function Yv = Y(obj, s)
            Yv = 2 * pi * (1 - cos(obj.a)) * sqrt(obj.f / obj.rho) * s^2.5; 
        end

        function zv = z(obj, s, omega)
            zv = (2/3) * omega * sqrt(obj.rho * obj.f()) * s^1.5;
        end

        function [B_A, Yeff, A, B] = backpropagate(obj, omega, Yeff_2)
            [Yeff, A, B] = deal(0, 0, 0);
            
            switch (obj.type)
                case 2
                    %Impedance Z 
                    Z = (obj.RW1 + obj.RW2 - ...
                        1i*omega*obj.RW1*obj.RW2*obj.Cwk)....
                        / (1-1i*omega*obj.RW2*obj.Cwk);
                
                    %Bessel functions at s1out
                    s1out = obj.s(obj.L);
                    [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] =...
                        besselfunctions(obj.a,s1out,omega,obj.rho,obj.beta); 
                    Y_s1out = obj.Y(s1out);
                    
                    B_A = -(J_13_s1out+1i*Y_s1out*J_43s1out*Z)/...
                        (Y_13s1out+1i*Y_s1out*Y_43s1out*Z);
                    
                    %Bessel functions at s1in
                    s1in = obj.s(0);
                    [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = ...
                        besselfunctions(obj.a,s1in,omega,obj.rho,obj.beta);
                    Y_s1in = (2*pi*(1-cos(obj.a)))*(fs1in/obj.rho)^0.5*s1in^2.5;
                    Yeff = -1i*Y_s1in*(J_43s1in+B_A*Y_43s1in)/(J_13_s1in+B_A*Y_13s1in);
                case 3
                    s4out = obj.s(obj.L);
                    [J_13s4out,Y_13s4out,J_43s4out,Y_43s4out,fs4out] = besselfunctions(obj.a,s4out,omega,obj.rho,obj.beta);    
                    Y_s4out = (2*pi*(1-cos(obj.a)))*(fs4out/obj.rho)^0.5*s4out^2.5;
                    B_A = -(1i*Y_s4out*J_43s4out+Yeff_2*J_13s4out)/(Yeff_2*Y_13s4out+1i*Y_s4out*Y_43s4out);
                    
                    s4in = obj.s(0);
                    [J_13s4in,Y_13s4in,J_43s4in,Y_43s4in,fs4in] = besselfunctions(obj.a,s4in,omega,obj.rho,obj.beta);      
                    Y_s4in = (2*pi*(1-cos(obj.a)))*(fs4in/obj.rho)^0.5*s4in^2.5;
                    Yeff = -1i*Y_s4in*(J_43s4in+B_A*Y_43s4in)/(J_13s4in+B_A*Y_13s4in);
                case 5
                    s4out = obj.s(obj.L);
                    [J_13s4out,Y_13s4out,J_43s4out,Y_43s4out,fs4out] = besselfunctions(obj.a,s4out,omega,obj.rho,obj.beta);    
                    Y_s4out = (2*pi*(1-cos(obj.a)))*(fs4out/obj.rho)^0.5*s4out^2.5;
                    B_A = -(1i*Y_s4out*J_43s4out+Yeff_2*J_13s4out)/(Yeff_2*Y_13s4out+1i*Y_s4out*Y_43s4out);

                    s1in = obj.s(0);
                    [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(obj.a,s1in,omega,obj.rho,obj.beta); 
                    Y_s1in = (2*pi*(1-cos(obj.a)))*(fs1in/obj.rho)^0.5*s1in^2.5;
                    A = -1/(1i*Y_s1in*(s1in^-0.5)*(J_43s1in+B_A*Y_43s1in));
                    B = B_A*A;
            end
        end

        function [] = forwardpropagate(obj, omega)
            
        end
    end
end

