classdef TaperedAortaModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        L = 24.137*10^-2;                                   %Table 2 (Flores 2016)
        R = 1.5*10^-2;                                      %Table 2 (Flores 2016)
        rho = 1060;                                         %Table 2 (Flores 2016)
        be = (1-0.5^2)/(400*1.2);
        RW1 = 1.8503*10^7;                                  %Table 2 (Flores 2016)
        RW2 = 1.0492*10^8;                                  %Table 2 (Flores 2016)
        Cwk = 1.0163*10^-8;                                 %Table 2 (Flores 2016)

        % External Data
        UTA_BC = load('PreviousCode/automeris/uper_thoracic_aorta/UTA_BC.csv');
        Inlet_pressure = load('PreviousCode/automeris/tapered_aorta/Inlet_Pressure.csv');
        midpoint_pressure = load('PreviousCode/automeris/tapered_aorta/midpoint_pressure.csv');
        outlet_pressure = load('PreviousCode/automeris/tapered_aorta/outlet_pressure.csv');
        midpoint_flow = load('PreviousCode/automeris/tapered_aorta/midpoint_flow.csv');
        outlet_flow = load('PreviousCode/automeris/tapered_aorta/outlet_flow.csv');
    end
    
    methods
        function [sol, toterr] = model(obj)
            L = obj.L;                                   %Table 2 (Flores 2016)
            R = obj.R;                                      %Table 2 (Flores 2016)
            rho = obj.rho;                                         %Table 2 (Flores 2016)
            RW1 = obj.RW1;                                  %Table 2 (Flores 2016)
            RW2 = obj.RW2;                                  %Table 2 (Flores 2016)
            Cwk = obj.Cwk;                                 %Table 2 (Flores 2016)
            be = obj.be;

            a = atan((R-10^-2)/L);

            x = obj.UTA_BC(1:end,1);
            y = obj.UTA_BC(1:end,2)*10^-6;
            N = length(y);
            T = obj.UTA_BC(end,1);
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
                  omega(ih+1) = ih*om;
                  
              end
              ih=ih+1;
                %% Step 1
                %Firstly the term B1/A1 
                
                %Impedance Z 
                Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);
              
                %Bessel functions at s1out
                s1out = (R-tan(a)*L)/sin(a);
                [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),rho,be); 
                Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho)^0.5*s1out^2.5;
                  
                
                B1_A1 = -(J_13_s1out+1i*Y_s1out*J_43s1out*Z)/(Y_13s1out+1i*Y_s1out*Y_43s1out*Z);
                
                s1in = (R-tan(a)*0)/sin(a);
                [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a,s1in,omega(ih),rho,be); 
                Y_s1in = (2*pi*(1-cos(a)))*(fs1in/rho)^0.5*s1in^2.5;
                A1tilda = -1/(1i*Y_s1in*(s1in^-0.5)*(J_43s1in+B1_A1*Y_43s1in));
                B1tilda = B1_A1*A1tilda;
                
                x1 = 0;
                s1 = (R-tan(a)*x1)/sin(a);
                [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),rho,be); 
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
                Q1(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
                P1(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
                
                x1 = L/2;
                s1 = (R-tan(a)*x1)/sin(a);
                [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),rho,be); 
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
                Q1mid(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
                P1mid(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
                
                x1 = L;
                s1 = (R-tan(a)*x1)/sin(a);
                [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),rho,be); 
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
                Q1outlet(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
                P1outlet(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
                
            end
            
            %% Inverse Fourier
            q1 = zeros(size(t)) + real(Q1(1)*F(1));
            p1 = zeros(size(t)) + real(P1(1)*F(1));
            q1mid = zeros(size(t)) + real(Q1mid(1)*F(1));
            p1mid = zeros(size(t)) + real(P1mid(1)*F(1));
            q1outlet = zeros(size(t)) + real(Q1outlet(1)*F(1));
            p1outlet = zeros(size(t)) + real(P1outlet(1)*F(1));
            
             for ih = 1:nh
                 
                    
                  q1 = q1 + real(conj(Q1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
                  p1 = p1 + real(conj(P1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
                    
                  q1mid = q1mid + real(conj(Q1mid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
                  p1mid = p1mid + real(conj(P1mid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
                  
                  q1outlet = q1outlet + real(conj(Q1outlet(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
                  p1outlet = p1outlet + real(conj(P1outlet(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
            
             end
            sol = [P1; Q1mid; P1mid; Q1outlet; P1outlet];
             
            
            %% Error calculations
            %Average
            yiinpressure = interp1q(obj.Inlet_pressure(1:end,1),obj.Inlet_pressure(1:end,2),xi);
            for i = 2:length(t)-2
                
                errinlet(i) = abs((yiinpressure(i) - p1(i))/(yiinpressure(i)))^2;
            
            end
            errorinletpressure = mean(errinlet);
            
            yimidpressure = interp1q(obj.midpoint_pressure(1:end,1),obj.midpoint_pressure(1:end,2),xi);
            for i = 2:length(t)-2
                
                errmid(i) = abs((yimidpressure(i) - p1mid(i))/(yimidpressure(i)))^2;
            
            end
            errormidpressure = mean(errmid);
            
            yimidflow = interp1q(obj.midpoint_flow(1:end,1),obj.midpoint_flow(1:end,2),xi);
            for i = 2:length(t)-3
                
                 errmid(i) = abs((yimidflow(i)*10^-6 - q1mid(i))/(max(yimidflow)*10^-6))^2;
            
            end
            errormidflow = mean(errmid);
            
            yioutletpressure = interp1q(obj.outlet_pressure(1:end,1),obj.outlet_pressure(1:end,2),xi);
            for i = 2:length(t)-2
                
                errmid(i) = abs((yioutletpressure(i) - p1outlet(i))/(yioutletpressure(i)))^2;
            
            end
            erroroutletpressure = mean(errmid);
            
            yioutletflow = interp1q(obj.outlet_flow(1:end,1),obj.outlet_flow(1:end,2),xi);
            for i = 2:length(t)-2
                
                 errmid(i) = abs((yioutletflow(i)*10^-6 - q1outlet(i))/(max(yioutletflow)*10^-6))^2;
            
            end
            erroroutletflow = mean(errmid);

            toterr = erroroutletflow + erroroutletpressure + errormidflow +...
                errormidpressure + errorinletpressure;
        end

        function err = globalerr(obj, xp, scaler, params, actsol)
            x = scaler.inv_transform(xp);
            for field = 1:length(params)
                obj.(params(field)) = x(field);
            end
            sol = obj.model();
            spe = (abs(sol - actsol) .^ 2) ./ abs(actsol) .^ 2 ;
            mspe = mean(spe, 2);
            err = sum(mspe);
        end

        function err = measurementerr(obj, xp, scaler, params)
            x = scaler.inv_transform(xp);
            for field = 1:length(params)
                obj.(params(field)) = x(field);
            end
            [~, err] = obj.model();
        end
    end
end

