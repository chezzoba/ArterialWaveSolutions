%% This is a function that implements the type III bifurcation
%This kind of bifurcation has one outlet B.C. (Windkesel Model).
%It requires an input (effective admitance) from the parent vessel of the 
%previously calculated bifurcation. This parent vessel is now a daughter
%vessel.
%It provides as output, information for the parent vessel.

function [B1_A1,B2_A2,Yeff_1] = type_III_bifurcation(Yeff_3,R1,R2,L1,L2,RW21,RW22,CWK2,beta1,beta2,rho1,rho2,omega,a1,a2)


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
    
    %% Parent Vessel --> 1
    %Calculation B1/A1
    s1out = (R1-tan(a1)*L1)/sin(a1);
    [J_13s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a1,s1out,omega,rho1,beta1);      
    Y_s1out = (2*pi*(1-cos(a1)))*(fs1out/rho1)^0.5*s1out^2.5;
    B1_A1 = -(1i*Y_s1out*J_43s1out+(Yeff_s2in+Yeff_3)*J_13s1out)/((Yeff_s2in+Yeff_3)*Y_13s1out+1i*Y_s1out*Y_43s1out);
    
    %Calculation of Yeff(s1in)
    s1in = (R1-tan(a1)*0)/sin(a1);
    [J_13s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a1,s1in,omega,rho1,beta1);      
    Y_s1in = (2*pi*(1-cos(a1)))*(fs1in/rho1)^0.5*s1in^2.5;
    Yeff_1 = -1i*Y_s1in*(J_43s1in+B1_A1*Y_43s1in)/(J_13s1in+B1_A1*Y_13s1in);
    
end