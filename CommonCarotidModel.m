classdef CommonCarotidModel
    %COMMONCAROTIDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a = 0.001;                                          %Very small angle
        %Properties
        L = 126*10^-3;                                      %Table 1 (Flores 2016)
        R = 3*10^-3;                                        %Table 1 (Flores 2016)
        rho = 1060;                                         %Table 1 (Flores 2016)  
        RW1 = 2.4875*10^8;                                  %Table 1 (Flores 2016)
        RW2 = 1.8697*10^9;                                  %Table 1 (Flores 2016)
        Cwk = 1.7529*10^-10;                                %Table 1 (Flores 2016)
        be = 0.003571428571429;
        % External Data
        CC_inlet_BC = load('PreviousCode/automeris/common_carotid_artery/CC_inlet_BC.csv');
        inlet_pressure = load('PreviousCode/automeris/common_carotid_artery/inlet_pressure.csv');
        mid_pressure = load('PreviousCode/automeris/common_carotid_artery/mid_pressure.csv');
        outlet_pressure = load('PreviousCode/automeris/common_carotid_artery/outlet_pressure.csv');
        CC_mid = load('PreviousCode/grabit/CC_mid.mat').CC_mid;
        CC_outlet = load('PreviousCode/grabit/CC_outlet.mat').CC_outlet;
    end
    
    methods
        function sol = model(obj)
            a = obj.a;
            %Properties
            L = obj.L;
            R = obj.R;
            rho = obj.rho;     
            RW1 = obj.RW1;
            RW2 = obj.RW2;
            Cwk = obj.Cwk;
            E = 5e5;
            v = 0.5;
            h = (1 - v^2) / (E * obj.be);

            %% Importing the data from the Flores plots
            CC_inlet_BC = obj.CC_inlet_BC;
            
            x = CC_inlet_BC(1:end,1);                            %Extracting the time data to a vector
            y = CC_inlet_BC(1:end,2)*10^-6;                      %Extracting the Q data to a vector
            N = length(y);
            T = CC_inlet_BC(end,1);                              %Extracting the period of th signal
            xi = (x(1):T/((N)):T)';                              %Creating an equispaced time vector to be used in the fft
            yi = interp1q(x,y,xi);                               %Interpolating the values of the Q to match the new time vector
            Qin = yi(1:end);                                     %Vector with the values of the BC
            t = xi(1:end);
            N = length(Qin);
            
            %% Fourier Transform
            F = fft(Qin(1:N))/N;
            
            % Define the fundamental harmonic omega and the number of harmonics
            om = 2*pi/T;
            nh = N/2;
            
            %% Frequency domain calculations
            for ih = 0:nh  
                 
              if ih == 0 
                  
                  omega(ih+1) = 10^-10;                               %Assining a very small value to ω to avoid singularity
                  
              else
                  omega(ih+1) = -ih*om;                               
                  
              end
              ih=ih+1;
            %% Backward Propagation 
                
                %Impedance Z 
                Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);
            
                %Bessel functions at s1out
                s1out = (R-tan(a)*L)/sin(a);
                [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),E,rho,v,h); 
                Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho)^0.5*s1out^2.5;
                
                B1_A1 = -(J_13_s1out+1i*Y_s1out*J_43s1out*Z)/(Y_13s1out+1i*Y_s1out*Y_43s1out*Z);
                
                s1in = (R-tan(a)*0)/sin(a);
                [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a,s1in,omega(ih),E,rho,v,h); 
                Y_s1in = (2*pi*(1-cos(a)))*(fs1in/rho)^0.5*s1in^2.5;
                A1tilda = -1/(1i*Y_s1in*(s1in^-0.5)*(J_43s1in+B1_A1*Y_43s1in));
                B1tilda = B1_A1*A1tilda;
              
            %% Forward Propagation     
                x1 = 0;
                s1 = (R-tan(a)*x1)/sin(a);
                [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E,rho,v,h); 
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
                Q1(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
                P1(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
                
                x1 = L/2;
                s1 = (R-tan(a)*x1)/sin(a);
                [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E,rho,v,h); 
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
                Q1mid(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
                P1mid(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
                
                x1 = L;
                s1 = (R-tan(a)*x1)/sin(a);
                [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E,rho,v,h); 
                Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
                Q1outlet(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
                P1outlet(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
            end
            sol = [P1; Q1mid; P1mid; Q1outlet; P1outlet];
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
    end
end

