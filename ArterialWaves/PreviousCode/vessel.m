function [Q,P,A,B] = vessel(Pout,L,R,a,omega,E,rho,v,h,B_A)

    xin = 0;
    s1in = (R-tan(a)*xin)/sin(a);
    [J_13sin,Y_13sin,J_43sin,Y_43sin,fsin] = besselfunctions(a,s1in,omega,E,rho,v,h); 
    A = Pout/((s1in^-0.5)*(J_13sin+B_A*Y_13sin));
    B = B_A*A;
    
    xout = L;
    sout = (R-tan(a)*xout)/sin(a);
    [J_13sout,Y_13out,J_43sout,Y_43sout,fsout] = besselfunctions(a,sout,omega,E,rho,v,h); 
    Y_sout = (2*pi*(1-cos(a)))*(fsout/rho)^0.5*sout^2.5;
    Q = -1i*Y_sout*(sout^-0.5)*(A*J_43sout+B*Y_43sout);
    P = ((sout^-0.5)*(A*J_13sout+B*Y_13out));
    
end