classdef BifurcationModel
    %BIFURCATIONMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Parent vessel (1)
        be1 = 0.75/(500*1.032);
        R1 = 0.89*10^-2;
        L1 = 8.6*10^-2;
 
        
        %Daughter vessels (2) and (3)
        be2 = 0.75/(700*0.72);
        R2 = 0.6125*10^-2;
        L2 = 8.5*10^-2;
        
        
        be3 = 0.75/(700*0.72);
        R3 = 0.6125*10^-2;
        L3 = 8.5*10^-2;
        
        %Impedence calculation data
        RW1 = 6.8123*10^7;
        RW2 = 3.1013*10^9;
        Cwk = 3.6664*10^-10;
        
        RW13 = 6.8123*10^7;
        RW23 = 3.1013*10^9;
        Cwk3 = 3.6664*10^-10;

        % Data
        ao_mid_flow = load('PreviousCode/grabit/ao_mid_flow.mat').ao_mid_flow;
        Inlet_ao_Pressure = load('PreviousCode/automeris/bifurcation/Inlet_ao_Pressure.csv');
        mid_ao_Pressure = load('PreviousCode/automeris/bifurcation/mid_ao_Pressure.csv');
        junction_Pressure = load('PreviousCode/automeris/bifurcation/junction_Pressure.csv');
        mid_il_Pressure = load('PreviousCode/automeris/bifurcation/mid_il_Pressure.csv');
        outlet_il_Pressure = load('PreviousCode/automeris/bifurcation/outlet_il_Pressure.csv');
        BC = load('PreviousCode/automeris/bifurcation/BC.csv');
        junction_flow = load('PreviousCode/automeris/bifurcation/junction_flow.csv');
        mid_il_flow = load('PreviousCode/automeris/bifurcation/mid_il_flow.csv');
        outlet_il_flow = load('PreviousCode/automeris/bifurcation/outlet_il_flow.csv');
    end
    
    methods
        function [sol, toterr] = model(obj, glob)
            rho1 = 1060;
            rho2 = 1060;
            rho3 = 1060;
            a = 0.001;
            a2 = 0.001;

            L1 = obj.L1;
            R1 = obj.R1;
            beta1 = obj.be1;

            L2 = obj.L2;
            R2 = obj.R2;
            beta2 = obj.be2;

            L3 = obj.L3;
            R3 = obj.R3;
            beta3 = obj.be3;

            RW1 = obj.RW1;
            RW2 = obj.RW2;
            Cwk = obj.Cwk;
            RW13 = obj.RW13;
            RW23 = obj.RW23;
            Cwk3 = obj.Cwk3;


            c1 = (1/(beta1*rho1*2*R1))^0.5;
            Y1 = (pi*R1^2)/(rho1*c1);
            c2 = (1/(beta2*rho2*2*R2))^0.5;
            Y2 = (pi*R2^2)/(rho2*c2);
            Y3 = Y2;
            
            x = obj.BC(1:end,1);
            y = obj.BC(1:end,2)*10^-6;
            N = length(y);
            T = obj.BC(end,1);
            xi = (x(1):T/((N)):T)';
            yi = interp1q(x,y,xi);
            Qin = yi(1:end);
            t = xi(1:end);
            N = length(Qin);
            
            %% Fourier Transform
            F = fft(Qin(1:N))/N;
            
            
            % Define the fundamental harmonic omega and the number of harmonics
            om = 2*pi/T;
            nh = N/2;
            
            %% Frequency domain calculations
            for ih = 0:nh   
            %% Backward Propagation       
                 
                if ih == 0
                    omega(ih+1) = 10^-10;
                
                else
                     omega(ih+1) = -ih*om;
                end
                
                ih=ih+1;
                
                %% Step 1
                %Firstly the terms B2/A2 and B3/A3 need to be calculated
                
                %Impedance Z (same for both daughter vessels)
                Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);
                Z3 = (RW13+RW23-1i*omega(ih)*RW13*RW23*Cwk3)/(1-1i*omega(ih)*RW23*Cwk3);
            
                %Bessel functions at s2out and s3out
                s2out = (R2-tan(a2)*L2)/sin(a2);
                s3out = (R3-tan(a2)*L3)/sin(a2);
                [J_13_s2out,Y_13s2out,J_43s2out,Y_43s2out,fs2out] = besselfunctions(a2,s2out,omega(ih),rho2,beta2); 
                [J_13_s3out,Y_13s3out,J_43s3out,Y_43s3out,fs3out] = besselfunctions(a2,s3out,omega(ih),rho3,beta3); 
                
                Y_s2out = (2*pi*(1-cos(a2)))*(fs2out/rho2)^0.5*s2out^2.5;
                Y_s3out = (2*pi*(1-cos(a2)))*(fs3out/rho2)^0.5*s3out^2.5;
                B2_A2 = -(J_13_s2out+1i*Y_s2out*J_43s2out*Z)/(Y_13s2out+1i*Y_s2out*Y_43s2out*Z);
                B3_A3 = -(J_13_s3out+1i*Y_s3out*J_43s3out*Z3)/(Y_13s3out+1i*Y_s3out*Y_43s3out*Z3);
                
                
                %% Step 2
                %Calculation of Yeff(s2in) and Yeff(s3in)
                s2in = R2/sin(a2);
                s3in = R3/sin(a2);
            
                [J_13_s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),rho2,beta2); 
                [J_13_s3in,Y_13s3in,J_43s3in,Y_43s3in,fs3in] = besselfunctions(a2,s3in,omega(ih),rho3,beta3); 
            
                Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
                Y_s3in = (2*pi*(1-cos(a2)))*(fs3in/rho3)^0.5*s3in^2.5;
            
                Yeff_s2in = -1i*Y_s2in*(J_43s2in+B2_A2*Y_43s2in)/(J_13_s2in+B2_A2*Y_13s2in);
                Yeff_s3in = -1i*Y_s3in*(J_43s3in+B3_A3*Y_43s3in)/(J_13_s3in+B3_A3*Y_13s3in);
            
                
                %% Step 3
                %Calculation of B1/A1
                s1out = (R1-tan(a)*L1)/sin(a);
                [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),rho1,beta1);
                Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho1)^0.5*s1out^2.5;
                B1_A1  = -((Yeff_s2in+Yeff_s3in)*J_13_s1out+1i*Y_s1out*J_43s1out)/(1i*Y_s1out*Y_43s1out+(Yeff_s2in+Yeff_s3in)*Y_13s1out);
                
                %% Calculation of flow rate at a specific point
                x0 = 0;    %Inlet
                s0 = (R1-tan(a)*x0)/sin(a);
                [J_13_s0,Y_13s0,J_43s0,Y_43s0,fs0] = besselfunctions(a,s0,omega(ih),rho1,beta1);
                Y_s0 = (2*pi*(1-cos(a)))*(fs0/rho1)^0.5*s0^2.5;
                
                x1 = 0;    %Inlet of parent vessel
                s1 = (R1-tan(a)*x1)/sin(a);
                [J_13_s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),rho1,beta1);
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho1)^0.5*s1^2.5; 
                Q1(ih) = (Y_s1*(s1^-0.5)*(J_43s1+B1_A1*Y_43s1))/(Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
                P1(ih) = -((s1^-0.5)*((J_13_s1)+B1_A1*Y_13s1))/(1i*Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
                
                x1 = L1/2;    %Middle of parent vessel
                s1 = (R1-tan(a)*x1)/sin(a);
                [J_13_s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),rho1,beta1);
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho1)^0.5*s1^2.5;
                Q1mid(ih) = (Y_s1*(s1^-0.5)*(J_43s1+B1_A1*Y_43s1))/(Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
                P1mid(ih) = -((s1^-0.5)*((J_13_s1)+B1_A1*Y_13s1))/(1i*Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
                
                x1 = L1;    %Junction
                s1out = (R1-tan(a)*x1)/sin(a);
                [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),rho1,beta1);
                Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho1)^0.5*s1out^2.5;
                A1tilda = -1/(1i*Y_s0*(s0^(-0.5)*(J_43s0+B1_A1*Y_43s0)));
                B1tilda = B1_A1*A1tilda;
                P1junc(ih) = A1tilda*(s1out^-0.5)*((J_13_s1out)+B1_A1*Y_13s1out);
                Q1junc(ih) = -1i*(Y_s1out*(s1out^-0.5))*(A1tilda*J_43s1out+B1tilda*Y_43s1out);
            
            %% Forward Propagation
                
               %% Step 1
                s2in = R2/sin(a2);
                [J_13s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),rho2,beta2);
                Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
                A2tilda =  Q1junc(ih)/(2*(-1i*Y_s2in*(s2in^(-0.5))*(J_43s2in+B2_A2*Y_43s2in)));
                B2tilda = B2_A2*A2tilda;
               
               %% Step 2
               x2in = 0;    %Inlet at the daughter vessel
               s2in = (R2-tan(a2)*x2in)/sin(a2);
               [J_13s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),rho2,beta2);
               Y_s2in = (2*pi*(1-cos(a2)))*((fs2in/rho2)^0.5)*s2in^2.5;
               Q2in(ih) = A2tilda*(-1i*Y_s2in*(s2in^(-0.5))*(J_43s2in+B2_A2*Y_43s2in));
               P2in(ih) = A2tilda*((s2in^(-0.5))*(J_13s2in+B2_A2*Y_13s2in));
               
               x2 = L2/2;    %Midlle of the daughter vessel
               s2 = (R2-tan(a2)*x2)/sin(a2);
               [J_13s2,Y_13s2,J_43s2,Y_43s2,fs2] = besselfunctions(a2,s2,omega(ih),rho2,beta2);
               Y_s2 = (2*pi*(1-cos(a2)))*((fs2/rho2)^0.5)*s2^2.5;
               Q2ilmid(ih) = A2tilda*(-1i*Y_s2*(s2^(-0.5))*(J_43s2+B2_A2*Y_43s2));  
               P2ilmid(ih) = A2tilda*((s2^(-0.5))*(J_13s2+B2_A2*Y_13s2));
               
               x2 = L2;     %Outlet of the daughter vessel
               s2 = (R2-tan(a2)*x2)/sin(a2);
               [J_13_s2,Y_13s2,J_43s2,Y_43s2,fs2] = besselfunctions(a2,s2,omega(ih),rho2,beta2);
               Y_s2 = (2*pi*(1-cos(a2)))*((fs2/rho2)^0.5)*s2^2.5;
               Q2ilout(ih) = A2tilda*(-1i*Y_s2*(s2^(-0.5))*(J_43s2+B2_A2*Y_43s2));
               P2ilout(ih) = Q2ilout(ih)*Z;
               
            end

            sol = [P1; Q1mid; P1mid; Q1junc; P1junc; Q2in;...
                P2in; Q2ilmid; P2ilmid; Q2ilout; P2ilout];
            toterr = 0;
            if (~glob)
                %% Inverse Fourier
                q1 = zeros(size(t)) + real(F(1));
                q1mid = zeros(size(t)) + real(Q1mid(1)*F(1));
                q1junc = zeros(size(t)) + real(Q1junc(1)*F(1));
                q2ilmid = zeros(size(t)) + real(Q2ilmid(1)*F(1));
                q2ilout = zeros(size(t)) + real(Q2ilout(1)*F(1));
                p1in = zeros(size(t)) + real(P1(1)*F(1));
                p1mid = zeros(size(t)) + real(P1mid(1)*F(1));
                p1junc = zeros(size(t)) + real(P1junc(1)*F(1));
                p2in = zeros(size(t)) + real(P2in(1)*F(1));
                p2mid = zeros(size(t)) + real(P2ilmid(1)*F(1));
                p2ilout = zeros(size(t)) + real(P2ilout(1)*F(1));
                
                
                 for ih = 1:nh
                     
                       q1 = q1 + real(Q1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       q1mid = q1mid + real(Q1mid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       q1junc = q1junc + real(Q1junc(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       q2ilmid = q2ilmid + real(Q2ilmid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       q2ilout = q2ilout + real(Q2ilout(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       p1in = p1in + real(P1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       p1mid = p1mid + real(P1mid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       p1junc = p1junc + real(P1junc(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       p2in = p2in + real(P2in(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       p2mid = p2mid + real(P2ilmid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                       p2ilout = p2ilout + real(P2ilout(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
                   
                 end
                
                %% Error calculations
                %Average
                yiinletpressure = interp1q(obj.Inlet_ao_Pressure(1:end,1),obj.Inlet_ao_Pressure(1:end,2),xi);
                for i = 2:length(t)-2
                    
                    errinletp(i) = abs((yiinletpressure(i) - p1in(i))/(yiinletpressure(i)))^2;
                
                end
                errorinletpressure = (mean(errinletp));
                
                
                yiaomidpressure = interp1q(obj.mid_ao_Pressure(1:end,1),obj.mid_ao_Pressure(1:end,2),xi);
                for i = 2:length(t)-2
                    
                    errmidp(i) = abs((yiaomidpressure(i) - p1mid(i))/(yiaomidpressure(i)))^2;
                
                end
                erroraomidpressure = (mean(errmidp));
                
                yiaomidflow = interp1q(obj.ao_mid_flow(1:end,1),obj.ao_mid_flow(1:end,2),xi);
                for i = 2:length(t)-2
                    
                     errmidf(i) = abs((yiaomidflow(i)*10^-6 - q1mid(i))/(max(yiaomidflow)*10^-6))^2;
                
                end
                erroraomidflow = (mean(errmidf));
                
                yijuncpressure = interp1q(obj.junction_Pressure(1:end,1),obj.junction_Pressure(1:end,2),xi);
                for i = 2:length(t)-2
                    
                    errjuncp(i) = abs((yijuncpressure(i) - p1junc(i))/(yijuncpressure(i)))^2;
                
                end
                errorjuncpressure = (mean(errjuncp));
                
                yijuncflow = interp1q(obj.junction_flow(1:end,1),obj.junction_flow(1:end,2),xi);
                for i = 1:length(t)-12
                    
                     errjuncf(i) = abs((yijuncflow(i)*10^-6 - q1junc(i))/(max(yijuncflow)*10^-6))^2;
                
                end
                errorjuncflow = (mean(errjuncf));
                
                yiilmidpressure = interp1q(obj.mid_il_Pressure(1:end,1),obj.mid_il_Pressure(1:end,2),xi);
                for i = 2:length(t)-2
                    
                    errillmidp(i) = abs((yiilmidpressure(i) - p2mid(i))/(yiilmidpressure(i)))^2;
                
                end
                errorilmidpressure = (mean(errillmidp));
                
                yiilmidflow = interp1q(obj.mid_il_flow(1:end,1),obj.mid_il_flow(1:end,2),xi);
                for i = 3:length(t)-3
                    
                     errilmidf(i) = abs((yiilmidflow(i)*10^-6 - q2ilmid(i))/(max(yiilmidflow)*10^-6))^2;
                
                end
                errorilmidflow = (mean(errilmidf));
                
                yiiloutpressure = interp1q(obj.outlet_il_Pressure(1:end,1),obj.outlet_il_Pressure(1:end,2),xi);
                for i = 2:length(t)-2
                    
                    erriloutp(i) = abs((yiiloutpressure(i) - p2ilout(i))/(yiiloutpressure(i)))^2;
                
                end
                erroriloutpressure = (mean(erriloutp));
                
                yiiloutflow = interp1q(obj.outlet_il_flow(1:end,1),obj.outlet_il_flow(1:end,2),xi);
                for i = 7:length(t)-10
                    
                     erriloutf(i) = abs((yiiloutflow(i)*10^-6 - q2ilout(i))/(max(yiiloutflow)*10^-6))^2;
                
                end
                erroriloutflow = (mean(erriloutf));
                
                toterr = sum([erroriloutflow, erroriloutpressure, errorilmidflow,...
                    errorilmidpressure, errorjuncflow, errorjuncpressure, ...
                    erroraomidflow, erroraomidpressure, errorinletpressure]);
            end
        end

        function err = globalerr(obj, xp, scaler, params)
            actsol = obj.model(true);
            x = scaler.inv_transform(xp);
            for field = 1:length(params)
                obj.(params(field)) = x(field);
            end
            sol = obj.model(true);
            spe = (abs(sol - actsol) .^ 2) ./ abs(actsol) .^ 2 ;
            mspe = mean(spe, 2);
            err = sum(mspe);
        end

        function err = measurementerr(obj, xp, scaler, params)
            x = scaler.inv_transform(xp);
            for field = 1:length(params)
                obj.(params(field)) = x(field);
            end
            [~, err] = obj.model(false);
        end

    end
end

