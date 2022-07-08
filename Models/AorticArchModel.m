classdef AorticArchModel
    %AORTICARCHMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rho = 1060;
        RW17 = 59149000;
        RW27 = 1.0174e+09;
        Cwk7 = 9.0686e-10;
        RW18 = 51918000;
        RW28 = 1060800000;
        Cwk8 = 8.6974e-10;
        RW19 = 191515000;
        RW29 = 5.22129e+09;
        Cwk9 = 1.767e-10;
        RW110 = 98820000;
        RW210 = 1.30183e+9;
        Cwk10 = 7.0871e-10;

        BC = load('PreviousCode/automeris/aorta/BC.csv');
        Inlet_Pressure = load('PreviousCode/automeris/aorta/Inlet_Pressure.csv');
        bc_outlet_11_flow = load('PreviousCode/automeris/aorta/bc_outlet_11_flow.csv');
        bc_outlet_11_Pressure = load('PreviousCode/automeris/aorta/bc_outlet_11_Pressure.csv');
        lcca_outlet_12_flow = load('PreviousCode/automeris/aorta/lcca_outlet_12_flow.csv');
        lcca_outlet_12_Pressure = load('PreviousCode/automeris/aorta/lcca_outlet_12_Pressure.csv');
        lsub_outlet_13_flow = load('PreviousCode/automeris/aorta/lsub_outlet_13_flow.csv');
        lsub_outlet_13_Pressure = load('PreviousCode/automeris/aorta/lsub_outlet_13_Pressure.csv');
    end
    
    methods

        function [sols, toterr] = model(obj)
            rho = obj.rho; % kg/m^3
            seg1 = 0.5; % L1
            WKP7 = [obj.RW17, obj.RW27, obj.Cwk7];
            
            
            ves1 = Vessel(0.0152, 0.0704*seg1, 0.0184750925619521, 0.001325687943664,...
                rho, [0, 0, 0]);
            ves1.type = 5;
            
            R2 = ves1.R - tan(ves1.a)*seg1*ves1.L;
            
            bi2_8_3 = Bifurcation([R2, 0.00635, 0.0139], [0.0704*(1-seg1), 0.0340, 0.008],...
                [0.00132568794366357, 0.00192990582059595, 0.00140439444384107],...
                rho, [0.0184750925619521, 0.001, 0.02499479361892],...
                obj.RW18, obj.RW28, obj.Cwk8);
            
            bi3_9_4 = Bifurcation([0.0139 0.0036 0.0137], [0.0080 0.0340 0.0090],...
                [0.001404394443841,0.002421354408802,0.001412397459944], rho,...
                [0.024994793618920,1e-03,0.022218565326719],...
                obj.RW19, obj.RW29, obj.Cwk9);
            
            bi4_10_5 = Bifurcation([0.0137,0.0048,0.0135], [0.009,0.034,0.064737],...
                [0.001412397459944,0.002158149171271,0.001388888888889], rho,...
                [0.022218565326719,1e-03,0.018534417520118],...
                obj.RW110, obj.RW210, obj.Cwk10);
            
            seg6 = 0.5;
            ves6 = Vessel(0.0123, seg6*0.152, 0.015788161735825, 0.001392773178531,...
                rho, [0, 0, 0]);
            
            R7 = ves6.R - tan(ves6.a)*seg6*ves6.L;
            
            ves7 = Vessel(R7, 0.152*(1-seg6), 0.015788161735825, 0.001392773178531,...
                rho, WKP7);
            ves7.type = 2;
            
            [omegas, F, t] = Vessel.ProcessBC(obj.BC);
            
            %% Backward Propagation
            
            ves7 = ves7.backpropagate(omegas);
            
            ves6 = ves6.backpropagate(omegas, ves7);
            
            ves5 = bi4_10_5.vessel(3).backpropagate(omegas, ves6);
            
            bi4_10_5 = bi4_10_5.backpropagate(omegas, ves5);
            
            bi3_9_4 = bi3_9_4.backpropagate(omegas, bi4_10_5);
            
            bi2_8_3 = bi2_8_3.backpropagate(omegas, bi3_9_4);
            
            ves1 = ves1.backpropagate(omegas, bi2_8_3);
            
            %% Forward Propagation
            
            %Inlet of vessel 1
            [Q1, P1] = ves1.forwardpropagate(omegas, ves1.s(0));
            
            %Outlet of vessel 1
            [Q1out, P1out] = ves1.forwardpropagate(omegas, ves1.s(ves1.L));
            
            ves2 = bi2_8_3.vessel(1);
            
            [Q2out, P2out] = ves2.forwardpropagate(omegas, ves2.s(ves2.L), P1out);
            
            [Q8out, P8out, Q3out, P3out] = bi2_8_3.forwardpropagate(omegas, P2out);
            
            [Q9out, P9out, Q4out, P4out] = bi3_9_4.forwardpropagate(omegas, P3out);
            
            [Q10out, P10out, Q5out, P5out] = bi4_10_5.forwardpropagate(omegas, P4out);
            
            [Q6out, P6out] = ves6.forwardpropagate(omegas, ves6.s(ves6.L), P5out);
            
            [Q7out, P7out] = ves7.forwardpropagate(omegas, ves7.s(ves7.L), P6out);
            
            
            %% Inverse Fourier Transform
            sols = [Q1; P1; Q8out; P8out; Q9out; P9out; Q10out; P10out];
            
            solt = InverseFourierTransform(t, omegas, sols, F);
            
            
            [q1, p1, q11, p11, q12, p12] = deal(solt(1, :), solt(2, :), solt(3, :), solt(4, :), solt(5, :), solt(6, :));
            
            [q13, p13] = deal(solt(7, :), solt(8, :));
            
            xi = t;
            %% Error calculations
            %Average
            
            yilsuboutletpressure = interp1q(obj.lsub_outlet_13_Pressure(1:end,1),obj.lsub_outlet_13_Pressure(1:end,2),xi);
            for i = 2:length(t)-3
                
                err13(i) = abs((yilsuboutletpressure(i) - p13(i))/(yilsuboutletpressure(i)))^2;
            
            end
            error13pressure = mean(err13);
            
            yilsuboutletflow = interp1q(obj.lsub_outlet_13_flow(1:end,1),obj.lsub_outlet_13_flow(1:end,2),xi);
            for i = 2:length(t)-4
                
                err13(i) = abs((yilsuboutletflow(i)*10^-6 - q13(i))/(max(yilsuboutletflow)*10^-6))^2;
            
            end
            error13flow = abs(mean(err13));
            
            yilccaoutletpressure = interp1q(obj.lcca_outlet_12_Pressure(1:end,1),obj.lcca_outlet_12_Pressure(1:end,2),xi);
            for i = 2:length(t)-3
                
                err12(i) = abs((yilccaoutletpressure(i) - p12(i))/(yilccaoutletpressure(i)))^2;
            
            end
            error12pressure = mean(err12);
            
            yilccaoutletflow = interp1q(obj.lcca_outlet_12_flow(1:end,1),obj.lcca_outlet_12_flow(1:end,2),xi);
            for i = 2:length(t)-2
                
                err12(i) = ((yilccaoutletflow(i)*10^-6 - q12(i))/(max(yilccaoutletflow)*10^-6))^2;
            
            end
            error12flow = abs(mean(err12));
            
            yibcoutletpressure = interp1q(obj.bc_outlet_11_Pressure(1:end,1),obj.bc_outlet_11_Pressure(1:end,2),xi);
            for i = 2:length(t)-3
                
                err11(i) = abs((yibcoutletpressure(i) - p11(i))/(yibcoutletpressure(i)))^2;
            
            end
            error11pressure = mean(err11);
            
            yibcoutletflow = interp1q(obj.bc_outlet_11_flow(1:end,1),obj.bc_outlet_11_flow(1:end,2),xi);
            for i = 3:length(t)-3
                
                err11(i) = ((yibcoutletflow(i)*10^-6 - q11(i))/(max(yibcoutletflow)*10^-6))^2;
            
            end
            error11flow = mean(err11);
            
            yiinletpressure = interp1q(obj.Inlet_Pressure(1:end,1),obj.Inlet_Pressure(1:end,2),xi);
            for i = 2:length(t)-3
                
                err1(i) = abs((yiinletpressure(i) - p1(i))/(yiinletpressure(i)))^2;
            
            end
            error1pressure = mean(err1);

            toterr = error1pressure + error11flow + error11pressure + error12flow...
                + error12pressure + error13flow + error13pressure;
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

