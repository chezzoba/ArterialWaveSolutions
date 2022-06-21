%% This is a function that implements the type II bifurcation
%This kind of bifurcation has two outlet B.C. (Windkesel).
%It does not require any input from other bifurcations
%It provides as output, information for the parent vessel.

function [B1_A1,B2_A2,B3_A3,Yeff_1] = type_II_bifurcation(R1,R2,R3,L1,L2,L3,RW21,RW22,CWK2,RW31,RW32,CWK3,beta1,beta2,beta3,rho1,rho2,rho3,omega,a1,a2,a3)
    
    %% Daughter vessel I --> 2
    %Calculation B2/A2
    Z2 = (RW21+RW22 - 1i*omega*RW21*RW22*CWK2)/(1-1i*omega*RW22*CWK2);
    s2out = (R2-tan(a2)*L2)/sin(a2);
    [J_13s2out,Y_13s2out,J_43s2out,Y_43s2out,fs2out] = besselfunctions(a2,s2out,omega,rho2,beta2);      
    Y_s2out = (2*pi*(1-cos(a2)))*(fs2out/rho2)^0.5*s2out^2.5;
    B2_A2 = -(J_13s2out+Z2*1i*Y_s2out*J_43s2out)/(1i*Z2*Y_s2out*Y_43s2out+Y_13s2out);
    
    %Calculation of Yeff(s2in)
    s2in = (R2-tan(a2)*0)/sin(a2);
    [J_13s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega,rho2,beta2);      
    Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
    Yeff_s2in = -1i*Y_s2in*(J_43s2in+B2_A2*Y_43s2in)/(J_13s2in+B2_A2*Y_13s2in);

    
    %% Daughter vessel II --> 3
    %Calculation B3/A3
    Z3 = (RW31+RW32 - 1i*omega*RW31*RW32*CWK3)/(1-1i*omega*RW32*CWK3);
    s3out = (R3-tan(a3)*L3)/sin(a3);
    [J_13s3out,Y_13s3out,J_43s3out,Y_43s3out,fs3out] = besselfunctions(a3,s3out,omega,rho3,beta3);      
    Y_s3out = (2*pi*(1-cos(a3)))*(fs3out/rho3)^0.5*s3out^2.5;
    B3_A3 = -(J_13s3out+Z3*1i*Y_s3out*J_43s3out)/(1i*Z3*Y_s3out*Y_43s3out+Y_13s3out);
    
    %Calculation of Yeff(s3in)
    s3in = (R3-tan(a3)*0)/sin(a3);
    [J_13s3in,Y_13s3in,J_43s3in,Y_43s3in,fs3in] = besselfunctions(a3,s3in,omega,rho3,beta3);      
    Y_s3in = (2*pi*(1-cos(a3)))*(fs3in/rho3)^0.5*s3in^2.5;
    Yeff_s3in = -1i*Y_s3in*(J_43s3in+B3_A3*Y_43s3in)/(J_13s3in+B3_A3*Y_13s3in);
    
    
    %% Parent Vessel --> 1
    %Calculation B1/A1
    s1out = (R1-tan(a1)*L1)/sin(a1);
    [J_13s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a1,s1out,omega,rho1,beta1);      
    Y_s1out = (2*pi*(1-cos(a1)))*(fs1out/rho1)^0.5*s1out^2.5;
    B1_A1 = -(1i*Y_s1out*J_43s1out+(Yeff_s2in+Yeff_s3in)*J_13s1out)/((Yeff_s2in+Yeff_s3in)*Y_13s1out+1i*Y_s1out*Y_43s1out);
    
    %Calculation of Yeff(s1in)
    s1in = (R1-tan(a1)*0)/sin(a1);
    [J_13s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a1,s1in,omega,rho1,beta1);      
    Y_s1in = (2*pi*(1-cos(a1)))*(fs1in/rho1)^0.5*s1in^2.5;
    Yeff_1 = -1i*Y_s1in*(J_43s1in+B1_A1*Y_43s1in)/(J_13s1in+B1_A1*Y_13s1in);
    
end