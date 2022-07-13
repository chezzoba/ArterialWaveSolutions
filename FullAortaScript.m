addpath("Algorithms/", "Models/", "VesselModels/");

aorta = FullAortaModel;
zfun = @(x, om) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
        ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]

measurements = 0;

switch (measurements)
    case 0
        [sols, ~, solt, t, omegas] = aorta.model();
        scaler = MinMaxScaler([1e7, 1e8, 0], [1e9,1e10,2e-9]);
        opt = NelderMeadSimplex;
        opt.lenX = 3;
        opt.x0Tol = 1e14;
        
        Z11 = sols(3, :) ./ sols(2, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z11).^2);
        WK11 = scaler.inv_transform(opt.fit(fun, 0))
        Perr11 = 100.*opt.relerr(WK11, [aorta.RW1(11), aorta.RW2(11), aorta.CWK(11)])

        Z12 = sols(5, :) ./ sols(4, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z12).^2);
        WK12 = scaler.inv_transform(opt.fit(fun, 0))
        Perr12 = 100.*opt.relerr(WK12, [aorta.RW1(12), aorta.RW2(12), aorta.CWK(12)])
        
        Z13 = sols(7, :) ./ sols(6, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z13).^2);
        WK13 = scaler.inv_transform(opt.fit(fun, 0))
        Perr13 = 100.*opt.relerr(WK13, [aorta.RW1(13), aorta.RW2(13), aorta.CWK(13)])

        Z15 = sols(9, :) ./ sols(8, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z15).^2);
        WK15 = scaler.inv_transform(opt.fit(fun, 0))
        Perr15 = 100.*opt.relerr(WK15, [aorta.RW1(15), aorta.RW2(15), aorta.CWK(15)])

        Z16 = sols(11, :) ./ sols(10, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z16).^2);
        WK16 = scaler.inv_transform(opt.fit(fun, 0))
        Perr16 = 100.*opt.relerr(WK16, [aorta.RW1(16), aorta.RW2(16), aorta.CWK(16)])

        Z18 = sols(13, :) ./ sols(12, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z18).^2);
        WK18 = scaler.inv_transform(opt.fit(fun, 0))
        Perr18 = 100.*opt.relerr(WK18, [aorta.RW1(18), aorta.RW2(18), aorta.CWK(18)])

        Z19 = sols(15, :) ./ sols(14, :);
        fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp), omegas) - Z19).^2);
        WK19 = scaler.inv_transform(opt.fit(fun, 0))
        Perr19 = 100.*opt.relerr(WK19, [aorta.RW1(19), aorta.RW2(19), aorta.CWK(19)])
        params = [];
        for i=1:20
            params = [params sprintf("be%d", i)];
        end

        scaler = MinMaxScaler(repmat(0.001, 1, 20), repmat(0.005, 1, 20));

        problem = OptimisationProblem(aorta, params, scaler);
        
        problem.optimiser.x0Tol = 13;
        %xpred1 = problem.fitsolution(0);
        %bePerr1 = 100 .* opt.relerr(problem.xact, scaler.inv_transform(xpred1));
        
        problem.optimiser = NelderMeadSimplex;
        problem.optimiser.MaxFunEvals = 20000;
        problem.optimiser.MaxIter = 20000;
        problem.optimiser.epochs = 1;
        %[xpred2, errP, erract, nguesses] = problem.fitsolution(xpred1);
        %bePerr2 = 100 .* opt.relerr(problem.xact, scaler.inv_transform(xpred2))
    case 1
        b = true;
end