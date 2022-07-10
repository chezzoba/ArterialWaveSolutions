function solt = InverseFourierTransform(ts, omegas, sols, F)
sols = conj(sols);   
solt = zeros(height(sols), length(F)) + real(sols(:, 1) .* F(1));
    
     for ih = 1:length(F)/2
         omt = repmat(omegas(ih+1) .* ts.', height(sols), 1);
         solt = solt + real(sols(:, ih+1) .* 2 .* F(ih+1) .* exp(1i*omt));
     end
end

