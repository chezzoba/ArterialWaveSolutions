function xpred = WindKessels(Pouts, Qouts)
    N = min([length(Pouts), length(Qouts)]) - 1;
    flow = interp1(Qouts(:, 1), Qouts(:, 2), Pouts(1:N, 1));
    
    Pout = fft(Pouts(1:N, 2)/N);
    Qout = fft(flow/N)*1e-6;
    om = -2*pi/(Qouts(N, 1)) * (0:N/2);
    om(1) = -1e-10;
    
    
    
    Z = Pout(1:floor(N/2+1)) ./ Qout(1:floor(N/2+1));
    
    zfun = @(x) (x(1) + x(2) - 1i .* om .* x(1) .* x(2) .* x(3))...
        ./ (1 - 1i .* om .* x(2) .* x(3)); % x = [RW1, RW2, Cwk]
    
    
    opt = ConOptimisation;
    opt.plt = false;
    scaler = MinMaxScaler([1e7, 1e8, 0], [1e9,1e10,2e-9]);
    
    fun = @(xp) mean(abs(zfun(scaler.inv_transform(xp)) - Z.').^2);
    opt.lenX = 3;
    opt.x0Tol = 1e30;
    xpred = scaler.inv_transform(opt.fit(fun, 0));
end

