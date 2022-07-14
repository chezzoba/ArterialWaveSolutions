addpath("Algorithms/", "Models/", "VesselModels/");
clc;
clear;

aorta = FullAortaModel;
zfun = @(x, om) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
        ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]
scaler = MinMaxScaler([1e7, 1e8, 0], [1e9,1e10,2e-9]);

measurements = 1;

switch (measurements)
    case 0
        [sols, ~, solt, t, omegas] = aorta.model();
        opt = NelderMeadSimplex;
        opt.lenX = 3;
        opt.x0Tol = 1e14;
        
%         Z11 = sols(3, :) ./ sols(2, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z11).^2);
%         WK11 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr11 = 100.*opt.relerr(WK11, [aorta.RW1(11), aorta.RW2(11), aorta.CWK(11)])
% 
%         Z12 = sols(5, :) ./ sols(4, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z12).^2);
%         WK12 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr12 = 100.*opt.relerr(WK12, [aorta.RW1(12), aorta.RW2(12), aorta.CWK(12)])
%         
%         Z13 = sols(7, :) ./ sols(6, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z13).^2);
%         WK13 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr13 = 100.*opt.relerr(WK13, [aorta.RW1(13), aorta.RW2(13), aorta.CWK(13)])
% 
%         Z15 = sols(9, :) ./ sols(8, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z15).^2);
%         WK15 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr15 = 100.*opt.relerr(WK15, [aorta.RW1(15), aorta.RW2(15), aorta.CWK(15)])
% 
%         Z16 = sols(11, :) ./ sols(10, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z16).^2);
%         WK16 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr16 = 100.*opt.relerr(WK16, [aorta.RW1(16), aorta.RW2(16), aorta.CWK(16)])
% 
%         Z18 = sols(13, :) ./ sols(12, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z18).^2);
%         WK18 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr18 = 100.*opt.relerr(WK18, [aorta.RW1(18), aorta.RW2(18), aorta.CWK(18)])
% 
%         Z19 = sols(15, :) ./ sols(14, :);
%         fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z19).^2);
%         WK19 = scaler.inv_transform(opt.fit(fun, 0))
%         Perr19 = 100.*opt.relerr(WK19, [aorta.RW1(19), aorta.RW2(19), aorta.CWK(19)])

        params = [];
        for i=1:20
            params = [params sprintf("be%d", i)];
        end
        
        params = [params ["RW114" "RW214" "CWK14"]];

        bescaler = MinMaxScaler([repmat(0.001, 1, 20) [1e7, 1e8, 0]],...
            [repmat(0.005, 1, 20) [1e9,1e10,2e-9]]);

        problem = OptimisationProblem(aorta, params, bescaler);
        
        problem.optimiser.x0Tol = 5;
        xpred1 = problem.fitsolution(0);
        bePerr1 = 100 .* opt.relerr(problem.xact, bescaler.inv_transform(xpred1));
        
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.MaxFunEvals = 40000;
        problem.optimiser.MaxIter = 40000;
        problem.optimiser.epochs = 1;
        [xpred2, errP, erract, nguesses] = problem.fitsolution(xpred1);
        bePerr2 = 100 .* opt.relerr(problem.xact, bescaler.inv_transform(xpred2))
    case 1
          opt = ConOptimisation;
          opt.lenX = 3;
          opt.x0Tol = 1e30;

        [q11, p11] = deal(aorta.bc_outlet_11_flow, aorta.bc_outlet_11_Pressure);
        N = length(p11) - 1;
        q11 = interp1(q11(:, 1), q11(:, 2), p11(1:N, 1));
        P11 = fft(p11(1:N, 2)/N);
        Q11 = fft(q11/N).*1e-6;
        om = -2*pi/(p11(N, 1)) * (0:N-1);
        Z11 = P11 ./ Q11;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z11.').^2);
        WK11 = scaler.inv_transform(opt.fit(fun, 0));
        Perr11 = 100.*opt.relerr(WK11, [aorta.RW1(11), aorta.RW2(11), aorta.CWK(11)])

        [q12, p12] = deal(aorta.lcca_outlet_12_flow, aorta.lcca_outlet_12_Pressure);
        N = length(p12) - 1;
        q12 = interp1(q12(:, 1), q12(:, 2), p12(2:N, 1));
        P12 = fft(p12(2:N, 2)/N);
        Q12 = fft(q12/N).*1e-6;
        om = -2*pi/(p12(N, 1)) * (0:N-2);
        Z12 = P12 ./ Q12;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z12.').^2);
        WK12 = scaler.inv_transform(opt.fit(fun, 0));
        Perr12 = 100.*opt.relerr(WK12, [aorta.RW1(12), aorta.RW2(12), aorta.CWK(12)])

        [q13, p13] = deal(aorta.lsub_outlet_13_flow, aorta.lsub_outlet_13_Pressure);
        N = length(p13) - 1;
        q13 = interp1(q13(:, 1), q13(:, 2), p13(1:N, 1));
        P13 = fft(p13(1:N, 2)/N);
        Q13 = fft(q13/N).*1e-6;
        om = -2*pi/(p13(N, 1)) * (0:N-1);
        Z13 = P13 ./ Q13;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z13.').^2);
        WK13 = scaler.inv_transform(opt.fit(fun, 0));
        Perr13 = 100.*opt.relerr(WK13, [aorta.RW1(13), aorta.RW2(13), aorta.CWK(13)])

        [q15, p15] = deal(aorta.sma_outlet_15_flow, aorta.sma_outlet_15_Pressure);
        N = length(p15) - 1;
        q15 = interp1(q15(:, 1), q15(:, 2), p15(2:N, 1));
        P15 = fft(p15(2:N, 2)/N);
        Q15 = fft(q15/N).*1e-6;
        om = -2*pi/(p15(N, 1)) * (0:N-2);
        Z15 = P15 ./ Q15;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z15.').^2);
        WK15 = scaler.inv_transform(opt.fit(fun, 0));
        Perr15 = 100.*opt.relerr(WK15, [aorta.RW1(15), aorta.RW2(15), aorta.CWK(15)])

        [q16, p16] = deal(aorta.rena_outlet_16_flow, aorta.rena_outlet_16_Pressure);
        N = length(p16) - 2;
        q16 = interp1(q16(:, 1), q16(:, 2), p16(1:N, 1));
        P16 = fft(p16(1:N, 2)/N);
        Q16 = fft(q16/N).*1e-6;
        om = -2*pi/(p16(N, 1)) * (0:N-1);
        Z16 = P16 ./ Q16;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z16.').^2);
        WK16 = scaler.inv_transform(opt.fit(fun, 0));
        Perr16 = 100.*opt.relerr(WK16, [aorta.RW1(16), aorta.RW2(16), aorta.CWK(16)])

        WK17 = WK16;

        [q18, p18] = deal(aorta.imma_outlet_18_flow, aorta.imma_outlet_18_Pressure);
        N = length(p18) - 1;
        q18 = interp1(q18(:, 1), q18(:, 2), p18(2:N, 1));
        P18 = fft(p18(2:N, 2)/N);
        Q18 = fft(q18/N).*1e-6;
        om = -2*pi/(p18(N, 1)) * (0:N-2);
        Z18 = P18 ./ Q18;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z18.').^2);
        WK18 = scaler.inv_transform(opt.fit(fun, 0));
        Perr18 = 100.*opt.relerr(WK18, [aorta.RW1(18), aorta.RW2(18), aorta.CWK(18)])

        [q19, p19] = deal(aorta.riliac_outlet_19_flow, aorta.riliac_outlet_19_Pressure);
        N = length(p19) - 1;
        q19 = interp1(q19(:, 1), q19(:, 2), p19(1:N, 1));
        P19 = fft(p19(1:N, 2)/N);
        Q19 = fft(q19/N).*1e-6;
        om = -2*pi/(p19(N, 1)) * (0:N-1);
        Z19 = P19 ./ Q19;
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), om) - Z19.').^2);
        WK19 = scaler.inv_transform(opt.fit(fun, 0));
        Perr19 = 100.*opt.relerr(WK19, [aorta.RW1(19), aorta.RW2(19), aorta.CWK(19)])

        WK20 = WK19;




end

