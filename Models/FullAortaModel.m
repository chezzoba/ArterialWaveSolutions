classdef FullAortaModel
    %FULLAORTAMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        be1 = 0.001325687943664;
        be2 = 0.001404394443841;
        be3 = 0.001412397459944;
        be4 = 0.001388888888889;
        be5 = 0.001392773178531;
        be6 = 0.001605713771886;
        be7 = 0.001624702408675;
        be8 = 0.001630675129943;
        be9 = 0.001614265805007;
        be10 = 0.001647214889066;
        be11 = 0.001929905820596;
        be12 = 0.002421354408802;
        be13 = 0.002158149171271;
        be14 = 0.002224647912390;
        be15 = 0.002382086707956;
        be16 = 0.002677500428400;
        be17 = 0.002677500428400;
        be18 = 0.003063224963241;
        be19 = 0.001973788094110;
        be20 = 0.001973788094110;
        RW1 = [0,0,0,0,0,0,0,0,0,0,5.1918,19.1515,9.882,11.7617,17.4352,34.1378,34.1378,74.0167,5.9149,5.9149]*(10^7);
        RW2 = [0,0,0,0,0,0,0,0,0,0,10.6080,52.2129,13.0183,7.5726,5.5097,5.3949,5.3949,46.2252,10.1737,10.1737]*(10^8);
        CWK = [0,0,0,0,0,0,0,0,0,0,8.6974,1.767,7.0871,12.1836,16.7453,17.1017,17.1017,1.9959,9.0686,9.0686]*(10^-10);
        RW114 = 117617000;
        RW214 = 757260000;
        CWK14 = 1.2184e-09;
        L = [7.0357,0.8,0.9,6.4737,15.2,1.8,0.7,0.7,4.3,4.3,3.4,3.4,3.4,3.2,6,3.2,3.2,5,8.5,8.5]*10^-2;
        Rin = [15.2,13.9,13.7,13.5,12.3,9.9,9.7,9.62,9.55,9.07,6.35,3.6,4.8,4.45,3.75,2.8,2.8,2,6,6]*10^-3;
        BC = load('PreviousCode/automeris/aorta/BC.csv');
        Inlet_Pressure = load('PreviousCode/automeris/aorta/Inlet_Pressure.csv');
        bc_outlet_11_flow = load('PreviousCode/automeris/aorta/bc_outlet_11_flow.csv');
        bc_outlet_11_Pressure = load('PreviousCode/automeris/aorta/bc_outlet_11_Pressure.csv');
        lcca_outlet_12_flow = load('PreviousCode/automeris/aorta/lcca_outlet_12_flow.csv');
        lcca_outlet_12_Pressure = load('PreviousCode/automeris/aorta/lcca_outlet_12_Pressure.csv');
        lsub_outlet_13_flow = load('PreviousCode/automeris/aorta/lsub_outlet_13_flow.csv');
        lsub_outlet_13_Pressure = load('PreviousCode/automeris/aorta/lsub_outlet_13_Pressure.csv');
        sma_outlet_15_flow = load('PreviousCode/automeris/aorta/sma_outlet_15_flow.csv');
        sma_outlet_15_Pressure = load('PreviousCode/automeris/aorta/sma_outlet_15_Pressure.csv');
        rena_outlet_16_flow = load('PreviousCode/automeris/aorta/rena_outlet_16_flow.csv');
        rena_outlet_16_Pressure = load('PreviousCode/automeris/aorta/rena_outlet_16_Pressure.csv');
        imma_outlet_18_flow = load('PreviousCode/automeris/aorta/imma_outlet_18_flow.csv');
        imma_outlet_18_Pressure = load('PreviousCode/automeris/aorta/imma_outlet_18_Pressure.csv');
        riliac_outlet_19_flow = load('PreviousCode/automeris/aorta/riliac_outlet_19_flow.csv');
        riliac_outlet_19_Pressure = load('PreviousCode/automeris/aorta/riliac_outlet_19_Pressure.csv');
    end
    
    methods
        
        function [sols, toterr, solt, t, omegas] = model(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            R=obj.Rin;
            L = obj.L;

            Rout = [R(2)*10^3,R(3)*10^3,R(4)*10^3,R(5)*10^3,R(6)*10^3,R(7)*10^3,R(8)*10^3,R(9)*10^3,...
            R(10)*10^3,8.6,R(11)*10^3,R(12)*10^3,R(13)*10^3,R(14)*10^3,R(15)*10^3,...
            R(16)*10^3,R(17)*10^3,R(18)*10^3,R(19)*10^3,R(20)*10^3]*10^-3;
            
            a(1:10) = atan((R(1:10)-Rout(1:10)) ./ L(1:10));
            a = [a(1:10),0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001];
            be = zeros(1, 20);
            for i=1:20
                be(i) = obj.(sprintf('be%d', i));
            end
            rho = 1060;
            
            [RW1, RW2, CWK] = deal(obj.RW1, obj.RW2, obj.CWK);

            bi10_19_20 = Bifurcation([R(10), R(19), R(20)], [L(10), L(19), L(20)],...
            [be(10), be(19), be(20)], rho, [a(10), a(19), a(20)],...
            [RW1(19), RW1(20)], [RW2(19), RW2(20)], [CWK(19), CWK(20)]);
        
            bi10_19_20.type = 2;
            
            bi9_18_10 = Bifurcation([R(9), R(18), R(10)], [L(9), L(18), L(10)],...
                [be(9), be(18), be(10)], rho, [a(9), a(18), a(10)],...
                RW1(18), RW2(18), CWK(18));
            
            bi8_17_9 = Bifurcation([R(8), R(17), R(9)], [L(8), L(17), L(9)],...
                [be(8), be(17), be(9)], rho, [a(8), a(17), a(9)],...
                RW1(17), RW2(17), CWK(17));
            
            bi7_16_8 = Bifurcation([R(7), R(16), R(8)], [L(7), L(16), L(8)],...
                [be(7), be(16), be(8)], rho, [a(7), a(16), a(8)],...
                RW1(16), RW2(16), CWK(16));
            
            bi6_15_7 = Bifurcation([R(6), R(15), R(7)], [L(6), L(15), L(7)],...
                [be(6), be(15), be(7)], rho, [a(6), a(15), a(7)],...
                RW1(15), RW2(15), CWK(15));
            
            bi5_14_6 = Bifurcation([R(5), R(14), R(6)], [L(5), L(14), L(6)],...
                [be(5), be(14), be(6)], rho, [a(5), a(14), a(6)],...
                obj.RW114, obj.RW214, obj.CWK14);
            
            bi3_13_4 = Bifurcation([R(3), R(13), R(4)], [L(3), L(13), L(4)],...
                [be(3), be(13), be(4)], rho, [a(3), a(13), a(4)],...
                RW1(13), RW2(13), CWK(13));
            
            bi2_12_3 = Bifurcation([R(2), R(12), R(3)], [L(2), L(12), L(3)],...
                [be(2), be(12), be(3)], rho, [a(2), a(12), a(3)],...
                RW1(12), RW2(12), CWK(12));
            
            bi1_11_2 = Bifurcation([R(1), R(11), R(2)], [L(1), L(11), L(2)],...
                [be(1), be(11), be(2)], rho, [a(1), a(11), a(2)],...
                RW1(11), RW2(11), CWK(11));
            
            bi1_11_2.type = 5;
            
            [omegas, F, t] = Vessel.ProcessBC(obj.BC);
               
            %% Backward Propagation  
            
            %10-19-20 --> Type II Bifurcation
            bi10_19_20 = bi10_19_20.backpropagate(omegas);
            
            %9-18-10 --> Type III Bifurcation
            bi9_18_10 = bi9_18_10.backpropagate(omegas, bi10_19_20);
            
            %8-17-9 --> Type III Bifurcation
            bi8_17_9 = bi8_17_9.backpropagate(omegas, bi9_18_10);
            
            %7-16-8 --> Type III Bifurcation
            bi7_16_8 = bi7_16_8.backpropagate(omegas, bi8_17_9);
            
            %6-15-7 --> Type III Bifurcation
            bi6_15_7 = bi6_15_7.backpropagate(omegas, bi7_16_8);
            
            %5-14-6 --> Type III Bifurcation
            bi5_14_6 = bi5_14_6.backpropagate(omegas, bi6_15_7);
            
            %4-5 --> Two connected vessels
            ves4 = bi3_13_4.vessel(3).backpropagate(omegas, bi5_14_6);
            
            %3-13-4 --> Type III Bifurcation  
            bi3_13_4 = bi3_13_4.backpropagate(omegas, ves4);
            
            %2-12-3 --> Type III Bifurcation    
            bi2_12_3 = bi2_12_3.backpropagate(omegas, bi3_13_4);
            
            %1-11-2 --> Type V Bifurcation 
            bi1_11_2 = bi1_11_2.backpropagate(omegas, bi2_12_3);
            
            %% Forward Propagation 
            %Inlet of vessel 1
            ves1 = bi1_11_2.vessel(1);
            
            [Q1, P1] = ves1.forwardpropagate(omegas, ves1.s(0));
            
            %Outlet of vessel 1
            [Q1out, P1out] = ves1.forwardpropagate(omegas, ves1.s(ves1.L));
            
            % Brachiocephalic and AO II
            [Q11out, P11out, Q2out, P2out] = bi1_11_2.forwardpropagate(omegas, P1out);
            
            %L com. carotid and AO III
            [Q12out, P12out, Q3out, P3out] = bi2_12_3.forwardpropagate(omegas, P2out);
            
            %Left subclavian and AO IV
            [Q13out, P13out, Q4out, P4out] = bi3_13_4.forwardpropagate(omegas, P3out);
            
            %AO V
            ves5 = bi5_14_6.vessel(1);
            [Q5out, P5out] = ves5.forwardpropagate(omegas, ves5.s(ves5.L), P4out);
            
            %AO VI
            [Q14out, P14out, Q6out, P6out] = bi5_14_6.forwardpropagate(omegas, P5out);
            
            %Sup. mesenteric and AO VII
            [Q15out, P15out, Q7out, P7out] = bi6_15_7.forwardpropagate(omegas, P6out);
            
            %Renal and AO VIII
            [Q16out, P16out, Q8out, P8out] = bi7_16_8.forwardpropagate(omegas, P7out);
            
            %AO IX and 17
            [Q17out, P17out, Q9out, P9out] = bi8_17_9.forwardpropagate(omegas, P8out);
            
            %Inf. mesenteric and AO X
            [Q18out, P18out, Q10out, P10out] = bi9_18_10.forwardpropagate(omegas, P9out);
            
            %R com. iliac
            [Q19out, P19out, Q20out, P20out] = bi10_19_20.forwardpropagate(omegas, P10out);
            
            
            % Inverse Fourier Transform
            sols = [P1; Q11out; P11out; Q12out; P12out; Q13out; P13out;
                Q15out; P15out; Q16out; P16out; Q18out; P18out; Q19out; P19out; Q1out];
            
            solt = Vessel.InverseFourierTransform(t, omegas, sols, F);
            
            xi = t;
            
            [p1, q11, p11, q12, p12] = deal(solt(1, :), solt(2, :), solt(3, :), solt(4, :), solt(5, :));
            [q13, p13, q15, p15, q16, p16] = deal(solt(6, :), solt(7, :), solt(8, :), solt(9, :), solt(10, :), solt(11, :));
            [q18, p18, q19, p19] = deal(solt(12, :), solt(13, :), solt(14, :), solt(15, :));
            
            %% Error calculations
            %Average
            yiliacoutletpressure = interp1q(obj.riliac_outlet_19_Pressure(1:end,1),obj.riliac_outlet_19_Pressure(1:end,2),xi);
            for i = 2:length(t)-2
                
                err19(i) = abs((yiliacoutletpressure(i) - p19(i))/(yiliacoutletpressure(i)))^2;
            
            end
            error19pressure = mean(err19);
            
            yiliacoutletflow = interp1q(obj.riliac_outlet_19_flow(1:end,1),obj.riliac_outlet_19_flow(1:end,2),xi);
            for i = 2:length(t)-2
                
                err19(i) = abs((yiliacoutletflow(i)*10^-6 - q19(i))/(max(yiliacoutletflow)*10^-6))^2;
            
            end
            error19flow = mean(err19);
            
            yiimmaoutletpressure = interp1q(obj.imma_outlet_18_Pressure(1:end,1),obj.imma_outlet_18_Pressure(1:end,2),xi);
            for i = 2:length(t)-2
                
                err18(i) = abs((yiimmaoutletpressure(i) - p18(i))/(yiimmaoutletpressure(i)))^2;
            
            end
            error18pressure = mean(err18);
            
            yiimmaoutletflow = interp1q(obj.imma_outlet_18_flow(1:end,1),obj.imma_outlet_18_flow(1:end,2),xi);
            for i = 2:length(t)-2
                
                err18(i) = abs((yiimmaoutletflow(i)*10^-6 - q18(i))/(max(yiimmaoutletflow)*10^-6))^2;
            
            end
            error18flow = mean(err18);
            
            yirenaoutletpressure = interp1q(obj.rena_outlet_16_Pressure(1:end,1),obj.rena_outlet_16_Pressure(1:end,2),xi);
            for i = 2:length(t)-2
                
                err16(i) = abs((yirenaoutletpressure(i) - p16(i))/(yirenaoutletpressure(i)))^2;
            
            end
            error16pressure = mean(err16);
            
            yirenaoutletflow = interp1q(obj.rena_outlet_16_flow(1:end,1),obj.rena_outlet_16_flow(1:end,2),xi);
            for i = 2:length(t)-5
                
                err16(i) = ((yirenaoutletflow(i)*10^-6 - q16(i))/(max(yirenaoutletflow)*10^-6))^2;
            
            end
            error16flow = mean(err16);
            
            yismaoutletpressure = interp1q(obj.sma_outlet_15_Pressure(1:end,1),obj.sma_outlet_15_Pressure(1:end,2),xi);
            for i = 2:length(t)-3
                
                err15(i) = abs((yismaoutletpressure(i) - p15(i))/(yismaoutletpressure(i)))^2;
            
            end
            error15pressure = mean(err15);
            
            yismaoutletflow = interp1q(obj.sma_outlet_15_flow(1:end,1),obj.sma_outlet_15_flow(1:end,2),xi);
            for i = 2:length(t)-5
                
                err15(i) = abs((yismaoutletflow(i)*10^-6 - q15(i))/(max(yismaoutletflow)*10^-6))^2;
            
            end
            error15flow = mean(err15);
            
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

            toterr = sum([error1pressure, error11flow, error11pressure, error12flow,...
                error12pressure, error13flow, error13pressure, error15flow, error15pressure,...
                error16flow, error16pressure, error18flow, error18pressure, error19flow, error19pressure]);
            
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
            [~, err] = obj.model(false);
        end

    end
end

