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

        data = readtable('Data/rsfs20170006supp1.xls', 'Sheet', sprintf('PD%d', 1));
        As = table2array(readtable('Data/rsfs20170006supp2.xls', 'Sheet', sprintf('PD%d', 1)));
        Qs = table2array(readtable('Data/rsfs20170006supp3.xls', 'Sheet', sprintf('PD%d', 1)));
    end
    
    methods

        function [sols, toterr] = model(obj)
            [As, Qs] = deal(obj.As(:, 2:end), obj.Qs(:, 2:end)*1e-6);
            rho = obj.rho; % kg/m^3
            Nx = height(Qs);

            avQ = mean(Qs, 2);
            Qdrops = find(abs((avQ(1:Nx-1) - avQ(2:Nx)) ./ avQ(1:Nx-1)) > 0.1);
            pos1 = Qdrops(1);
            pos3 = Qdrops(end);
            
            poses = round([1, pos1/2, pos1, (pos1 + pos3) / 2,...
                pos3, (Nx + 2*pos3)/3, (2*Nx + pos3)/3, Nx]);
            
            xs = obj.data.AnalysisPlane(poses) * 1e-3;
            Asysp = obj.data.MaximumAreaMeasured(poses) * 1e-6;
            Adiasp = obj.data.MinimumAreaMeasured(poses) * 1e-6;
            Psysp = obj.data.MaximumPressure(poses) * 133.32237;
            Pdiasp = obj.data.MinimumPressure(poses) * 133.32237;
            
            Ls = (xs(2:end) - xs(1:length(xs)-1));
            Rs = obj.data.EndDiastolicRadius(poses) * 1e-3;
            betas = (1 - sqrt(Adiasp ./ Asysp)) ./ (Rs .* (Psysp - Pdiasp));
            as = atan((Rs(1:length(Rs)-1) - Rs(2:end)) ./ Ls);
            
            parents = [2, 3, 4];
            chila = [3, 4, 5];
            Rpparents = sqrt(2 .* rho ./ betas(parents)) ./ (pi .* Rs(parents) .^ 2.5);
            Rpchila = sqrt(2 .* rho ./ betas(chila)) ./ (pi .* Rs(chila) .^ 2.5);
            Rpsa = (Rpparents .* Rpchila) ./ (Rpparents - Rpchila);
            betasa = 2 * rho ./ (Rpsa .* pi .* Rs(chila) .^ 2.5) .^ 2;
            Lsa = 2e-2;
            
            
            %% Importing the data from the Flores plots
            
            WKP7 = [obj.RW17, obj.RW27, obj.Cwk7];
            
            ves1 = Vessel(Rs(1), Ls(1), as(1), betas(1),...
                rho, [0, 0, 0]);
            ves1.type = 5;
            
            bi2_8_3 = Bifurcation([Rs(2), Rs(3), Rs(3)], [Ls(2), Lsa, Ls(3)],...
                [betas(2), betasa(1), betas(3)],...
                rho, [as(1), 0.001, as(3)],...
                obj.RW18, obj.RW28, obj.Cwk8);
            
            bi3_9_4 = Bifurcation([Rs(3) Rs(4) Rs(4)], [Ls(3) Lsa Ls(4)],...
                [betas(3),betasa(2),betas(4)], rho,...
                [as(3),1e-03,as(4)],...
                obj.RW19, obj.RW29, obj.Cwk9);
            
            bi4_10_5 = Bifurcation([Rs(4),Rs(5),Rs(5)], [Ls(4),Lsa,Ls(5)],...
                [betas(4),betasa(3),betas(5)], rho,...
                [as(4),1e-03,as(5)],...
                obj.RW110, obj.RW210, obj.Cwk10);
            
            ves6 = Vessel(Rs(6), Ls(6), as(6), betas(6), rho, [0, 0, 0]);
            
            ves7 = Vessel(Rs(7), Ls(7), as(7), betas(7), rho, WKP7);
            ves7.type = 2;

            xi = (0:24) / 25;
            BC = [xi; Qs(1, :)].';
            [omegas, F, t] = Vessel.ProcessBC(BC);
            
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
            sols = [Q1; P1; Q8out; P8out; Q9out; P9out; Q10out; P10out; Q7out; P7out];
            
            solt = Vessel.InverseFourierTransform(t, omegas, sols, F);
            
            
            [q1, p1, q8, p8, q9, p9] = deal(solt(1, :), solt(2, :), solt(3, :), solt(4, :), solt(5, :), solt(6, :));
            
            [q10, p10, q7, p7] = deal(solt(7, :), solt(8, :), solt(9, :), solt(10, :));
            q7meas = interp1(xi, Qs(end, :), t);
            toterr = mean(mean((q7meas - q7).^2));
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

