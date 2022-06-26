classdef Vessel
    %SINGLEVESSEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = 3; % 1, 2, 3 or 5
        R;
        L;
        a;
        beta;
        rho;
        RW1;
        RW2;
        Cwk;
        B1_A1 = 0;
        A1 = 0;
        oms = 0;
        Yeff = 0;
    end
    
    methods
        function obj = Vessel(R, L, a, be, rho, WKP)
            %SINGLEVESSEL Construct an instance of this class
            %   Detailed explanation goes here
            [obj.R, obj.L, obj.a, obj.beta, obj.rho] = ...
                deal(R, L, a, be, rho);
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
            Yv = 2 * pi * (1 - cos(obj.a)) * sqrt(obj.f() / obj.rho) * s^2.5; 
        end

        function zv = z(obj, s)
            zv = (2/3) * obj.oms * sqrt(obj.rho * obj.f()) * s^1.5;
        end

        function obj = backpropagate(obj, inp)

            switch (obj.type)
                case 1
                    obj.oms = inp;
                    %Impedance Z 
                    Z = (obj.RW1 + obj.RW2 - ...
                        1i .* obj.oms .* obj.RW1 .* obj.RW2 .* obj.Cwk)....
                        ./ (1-1i .* obj.oms .*obj.RW2 .* obj.Cwk);
                
                    %Bessel functions at s1out
                    s1out = obj.s(obj.L);
                    [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] =...
                        besselfunctions(obj.a,s1out,obj.oms,obj.rho,obj.beta); 
                    Y_s1out = obj.Y(s1out);
                    
                    obj.B1_A1 = -(J_13_s1out+1i .* Y_s1out .* J_43s1out .* Z) ./...
                        (Y_13s1out+1i .* Y_s1out .* Y_43s1out .* Z);

                    %Bessel functions at s1in
                    s1in = obj.s(0);
                    [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(obj.a,s1in,obj.oms,obj.rho,obj.beta); 
                    Y_s1in = obj.Y(s1in);
                    obj.A1 = -1 ./ (1i .* Y_s1in .* (s1in .^ -0.5) .* (J_43s1in+obj.B1_A1 .* Y_43s1in));
                case 2
                    %Impedance Z 
                    obj.oms = inp;
                    Z = (obj.RW1 + obj.RW2 - ...
                        1i*obj.oms*obj.RW1*obj.RW2*obj.Cwk)....
                        / (1-1i*obj.oms*obj.RW2*obj.Cwk);
                
                    %Bessel functions at s1out
                    s1out = obj.s(obj.L);
                    [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] =...
                        besselfunctions(obj.a,s1out,obj.oms,obj.rho,obj.beta); 
                    Y_s1out = obj.Y(s1out);
                    
                    obj.B1_A1 = -(J_13_s1out+1i*Y_s1out*J_43s1out*Z)/...
                        (Y_13s1out+1i*Y_s1out*Y_43s1out*Z);
                    
                    %Bessel functions at s1in
                    s1in = obj.s(0);
                    [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = ...
                        besselfunctions(obj.a,s1in,obj.oms,obj.rho,obj.beta);
                    Y_s1in = (2*pi*(1-cos(obj.a)))*(fs1in/obj.rho)^0.5*s1in^2.5;
                    obj.Yeff = -1i*Y_s1in*(J_43s1in+obj.B1_A1*Y_43s1in)/(J_13_s1in+obj.B1_A1*Y_13s1in);
                case 3
                    Yeff_2 = inp.Yeff;
                    obj.oms = inp.oms;
                    s4out = obj.s(obj.L);
                    [J_13s4out,Y_13s4out,J_43s4out,Y_43s4out,fs4out] = besselfunctions(obj.a,s4out,obj.oms,obj.rho,obj.beta);    
                    Y_s4out = (2*pi*(1-cos(obj.a)))*(fs4out/obj.rho)^0.5*s4out^2.5;
                    obj.B1_A1 = -(1i.*Y_s4out.*J_43s4out+Yeff_2.*J_13s4out)./(Yeff_2.*Y_13s4out+1i.*Y_s4out.*Y_43s4out);
                    
                    s4in = obj.s(0);
                    [J_13s4in,Y_13s4in,J_43s4in,Y_43s4in,fs4in] = besselfunctions(obj.a,s4in,obj.oms,obj.rho,obj.beta);      
                    Y_s4in = (2*pi*(1-cos(obj.a)))*(fs4in/obj.rho)^0.5*s4in^2.5;
                    obj.Yeff = -1i.*Y_s4in.*(J_43s4in+obj.B1_A1.*Y_43s4in)./(J_13s4in+obj.B1_A1.*Y_13s4in);
                case 5
                    Yeff_2 = inp.Yeff;
                    obj.oms = inp.oms;
                    s4out = obj.s(obj.L);
                    [J_13s4out,Y_13s4out,J_43s4out,Y_43s4out,fs4out] = besselfunctions(obj.a,s4out,obj.oms,obj.rho,obj.beta);    
                    Y_s4out = (2*pi*(1-cos(obj.a)))*(fs4out/obj.rho)^0.5*s4out^2.5;
                    obj.B1_A1 = -(1i.*Y_s4out.*J_43s4out+Yeff_2.*J_13s4out)./(Yeff_2.*Y_13s4out+1i.*Y_s4out.*Y_43s4out);

                    s1in = obj.s(0);
                    [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(obj.a,s1in,obj.oms,obj.rho,obj.beta); 
                    Y_s1in = (2*pi*(1-cos(obj.a)))*(fs1in/obj.rho)^0.5*s1in^2.5;
                    obj.A1 = -1./(1i.*Y_s1in.*(s1in^-0.5).*(J_43s1in+obj.B1_A1.*Y_43s1in));
            end
        end

        function [Q, P, A] = forwardpropagate(obj, s, P0outi)
            switch (obj.type)
                case {2, 3}
                    [Q, P, A] = vesselforward(P0outi,obj.L,obj.R,obj.a,obj.oms,obj.rho,obj.beta,obj.B1_A1);
                case {1, 5}
                    B = obj.B1_A1 .* obj.A1;
                    A = obj.A1;
                    [J_13s1,Y_13s1,J_43s1,Y_43s1,~] = besselfunctions(obj.a,s,obj.oms,obj.rho,obj.beta); 
                    Y_s1 = obj.Y(s);
                    Q = -(1i.*Y_s1.*(s^-0.5).*(A.*J_43s1+B.*Y_43s1));
                    P = ((s^-0.5).*(A.*J_13s1+B.*Y_13s1));
            end
        end
    end
end

