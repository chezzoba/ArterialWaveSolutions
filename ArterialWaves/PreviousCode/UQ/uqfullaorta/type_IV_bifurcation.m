%% This is a function that implements the type IV bifurcation
%This kind of bifurcation has no B.C.
%It requires two inputs (effective admitances) from the parent vessels of the 
%previously calculated bifurcations. These parent vessels are now the two
%daughter vessels.
%It provides as output, information for the parent vessel.

function [B1_A1,Yeff_1] = type_IV_bifurcation(Yeff_2,Yeff_3,R1,L1,E1,h1,v1,rho1,omega,a)


    %% Parent Vessel --> 1
    %Calculation B1/A1
    s1out = (R1-tan(a)*L1)/sin(a);
    [J_13s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega,E1,rho1,v1,h1);      
    Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho1)^0.5*s1out^2.5;
    B1_A1 = -(1i*Y_s1out*J_43s1out+(Yeff_2+Yeff_3)*J_13s1out)/((Yeff_2+Yeff_3)*Y_13s1out+1i*Y_s1out*Y_43s1out);
    
    %Calculation of Yeff(s1in)
    s1in = (R1-tan(a)*0)/sin(a);
    [J_13s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a,s1in,omega,E1,rho1,v1,h1);      
    Y_s1in = (2*pi*(1-cos(a)))*(fs1in/rho1)^0.5*s1in^2.5;
    Yeff_1 = -1i*Y_s1in*(J_43s1in+B1_A1*Y_43s1in)/(J_13s1in+B1_A1*Y_13s1in);
    
end