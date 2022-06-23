%% This is a function that implements the type IV bifurcation
%This kind of bifurcation has one inlet B.C and one outlet BC.
%It requires one input (effective admitance) from the parent vessel of the 
%previously calculated bifurcation. This parent vessel is now one
%daughter vessel.
%It provides as output, information for the parent vessel (A1,B1)and
%information about the daughter vessel (B2_A2)

function [A1,B1,B2_A2] = type_V_bifurcation(Yeff_3,R1,R2,L1,L2,RW21,RW22,CWK2,beta1,beta2,rho1,rho2,omega,a,a2)
    
%% Daughter vessel I --> 2
    %Calculation B2/A2
    Z2 = (RW21+RW22 - 1i.*omega.*RW21.*RW22.*CWK2)./(1-1i.*omega.*RW22.*CWK2);
    s2out = (R2-tan(a2)*L2)/sin(a2);
    [J_13s2out,Y_13s2out,J_43s2out,Y_43s2out,fs2out] = besselfunctions(a2,s2out,omega,rho2,beta2);      
    Y_s2out = (2*pi*(1-cos(a2)))*(fs2out/rho2)^0.5*s2out^2.5;
    B2_A2 = -(J_13s2out+Z2.*1i.*Y_s2out.*J_43s2out)./(1i.*Z2.*Y_s2out.*Y_43s2out+Y_13s2out);

    %Calculation of Yeff(s2in)
    s2in = (R2-tan(a2)*0)/sin(a2);
    [J_13s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega,rho2,beta2);      
    Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
    Yeff_2 = -1i.*Y_s2in.*(J_43s2in+B2_A2.*Y_43s2in)./(J_13s2in+B2_A2.*Y_13s2in);
    
    %% Parent Vessel --> 1
    %Calculation B1/A1
    s1out = (R1-tan(a)*L1)/sin(a);
    [J_13s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega,rho1,beta1);      
    Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho1)^0.5*s1out^2.5;
    B1_A1 = -(1i.*Y_s1out.*J_43s1out+(Yeff_2+Yeff_3).*J_13s1out)./((Yeff_2+Yeff_3).*Y_13s1out+1i.*Y_s1out.*Y_43s1out);

    
    %Calculation A1
    s1in = (R1-tan(a)*0)/sin(a);
    [J_13s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a,s1in,omega,rho1,beta1);      
    Y_s1in = (2*pi*(1-cos(a)))*(fs1in/rho1)^0.5*s1in^2.5;
    A1 = -1./(1i.*Y_s1in.*(s1in^-0.5).*(J_43s1in+B1_A1.*Y_43s1in));
    
    %Calculation B1
    B1 = B1_A1.*A1;

end