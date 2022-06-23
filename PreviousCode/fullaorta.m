%% Programme for implementation of the full aorta model
clear;
clc;
%% Flores et all (2016) and Nan Xiao et al (2013), Table 5 and Table IV, Parameters of the full-aorta model
L = [7.0357,0.8,0.9,6.4737,15.2,1.8,0.7,0.7,4.3,4.3,3.4,3.4,3.4,3.2,6,3.2,3.2,5,8.5,8.5]*10^-2;
Rin = [15.2,13.9,13.7,13.5,12.3,9.9,9.7,9.62,9.55,9.07,6.35,3.6,4.8,4.45,3.75,2.8,2.8,2,6,6]*10^-3;
Rout = [Rin(2)*10^3,Rin(3)*10^3,Rin(4)*10^3,Rin(5)*10^3,Rin(6)*10^3,Rin(7)*10^3,Rin(8)*10^3,Rin(9)*10^3,Rin(10)*10^3,8.6,Rin(11)*10^3,Rin(12)*10^3,Rin(13)*10^3,Rin(14)*10^3,Rin(15)*10^3,Rin(16)*10^3,Rin(17)*10^3,Rin(18)*10^3,Rin(19)*10^3,Rin(20)*10^3]*10^-3;
h = 0.1*Rin; %h was chosen to be 10% of Rin
R=Rin;
E = [372.2,384.2,387.6,400,437.8,471.8,475.9,478.1,486.5,502,612,860.4,724,757.6,839.6,1000.4,1000.4,1224.2,633.3,633.3]*10^3;
RW1 = [0,0,0,0,0,0,0,0,0,0,5.1918,19.1515,9.882,11.7617,17.4352,34.1378,34.1378,74.0167,5.9149,5.9149]*(10^7);
RW2 = [0,0,0,0,0,0,0,0,0,0,10.6080,52.2129,13.0183,7.5726,5.5097,5.3949,5.3949,46.2252,10.1737,10.1737]*(10^8);
CWK = [0,0,0,0,0,0,0,0,0,0,8.6974,1.767,7.0871,12.1836,16.7453,17.1017,17.1017,1.9959,9.0686,9.0686]*(10^-10);
v = 0.5;
rho = 1060; 
a(1:10) = atan((Rin(1:10)-Rout(1:10))./L(1:10));
a = [a(1:10),0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001];

%% Importing the data from the Flores plots
load('automeris/aorta/BC.csv')
load('automeris/aorta/Inlet_Pressure.csv')
load('automeris/aorta/bc_outlet_11_flow.csv')
load('automeris/aorta/bc_outlet_11_Pressure.csv')
load('automeris/aorta/lcca_outlet_12_flow.csv')
load('automeris/aorta/lcca_outlet_12_Pressure.csv')
load('automeris/aorta/lsub_outlet_13_flow.csv')
load('automeris/aorta/lsub_outlet_13_Pressure.csv')
load('automeris/aorta/sma_outlet_15_flow.csv')
load('automeris/aorta/sma_outlet_15_Pressure.csv')
load('automeris/aorta/rena_outlet_16_flow.csv')
load('automeris/aorta/rena_outlet_16_Pressure.csv')
load('automeris/aorta/imma_outlet_18_flow.csv')
load('automeris/aorta/imma_outlet_18_Pressure.csv')
load('automeris/aorta/riliac_outlet_19_flow.csv')
load('automeris/aorta/riliac_outlet_19_Pressure.csv')

tic
x = BC(1:end,1);
y = BC(1:end,2)*10^-6;
N = length(y);
T = BC(end,1);
xi = (x(1):T/((N)):T)';
yi = interp1q(x,y,xi);
Qin = yi(1:end);
t = xi(1:end);
N = length(Qin);

%% Fourier Transform
F = fft(Qin(1:N))/N;


% Define the fundamental harmonic omega and the number of harmonics
om = 2*pi/T;
nh = N/2;

% Memory allocation
omega = zeros(1,nh);
Q1 = zeros(1,nh);     P1 = zeros(1,nh);  
Q1out = zeros(1,nh);  P1out = zeros(1,nh);
Q2out = zeros(1,nh);  P2out = zeros(1,nh);  
Q3out = zeros(1,nh);  P3out = zeros(1,nh);
Q4out = zeros(1,nh);  P4out = zeros(1,nh);
Q5out = zeros(1,nh);  P5out = zeros(1,nh);
Q6out = zeros(1,nh);  P6out = zeros(1,nh);
Q7out = zeros(1,nh);  P7out = zeros(1,nh);
Q8out = zeros(1,nh);  P8out = zeros(1,nh);
Q9out = zeros(1,nh);  P9out = zeros(1,nh);
Q10out = zeros(1,nh); P10out = zeros(1,nh);
Q11out = zeros(1,nh); P11out = zeros(1,nh);
Q12out = zeros(1,nh); P12out = zeros(1,nh);
Q13out = zeros(1,nh); P13out = zeros(1,nh);
Q15out = zeros(1,nh); P15out = zeros(1,nh);
Q16out = zeros(1,nh); P16out = zeros(1,nh);
Q18out = zeros(1,nh); P18out = zeros(1,nh);
Q19out = zeros(1,nh); P19out = zeros(1,nh);

%% Calculations in the freequency domain
for ih = 0:nh     
    
    %% Freequency
    
    if ih == 0
         omega(ih+1) = 10^-10;
    
    else
         omega(ih+1) = ih*om;
    end
    
    ih=ih+1;
   
    %% Backward Propagation  

    %10-19-20 --> Type II Bifurcation
    [B10_A10,B19_A19,B20_A20,Yeff_10] = type_II_bifurcation(R(10),R(19),R(20),L(10),L(19),L(20),RW1(19),RW2(19),CWK(19),RW1(20),RW2(20),CWK(20),E(10),E(19),E(20),h(10),h(19),h(20),v,v,v,rho,rho,rho,omega(ih),a(10),a(19),a(20));
     
    %9-10-18 --> Type III Bifurcation
    [B9_A9,B18_A18,Yeff_9] = type_III_bifurcation(Yeff_10,R(9),R(18),L(9),L(18),RW1(18),RW2(18),CWK(18),E(9),E(18),h(9),h(18),v,v,rho,rho,omega(ih),a(9),a(18));
    
    %8-9-17 --> Type III Bifurcation
    [B8_A8,B17_A17,Yeff_8] = type_III_bifurcation(Yeff_9,R(8),R(17),L(8),L(17),RW1(17),RW2(17),CWK(17),E(8),E(17),h(8),h(17),v,v,rho,rho,omega(ih),a(8),a(17));
     
    %7-8-16 --> Type III Bifurcation
    [B7_A7,B16_A16,Yeff_7] = type_III_bifurcation(Yeff_8,R(7),R(16),L(7),L(16),RW1(16),RW2(16),CWK(16),E(7),E(16),h(7),h(16),v,v,rho,rho,omega(ih),a(7),a(16));
    
    %6-7-15 --> Type III Bifurcation
    [B6_A6,B15_A15,Yeff_6] = type_III_bifurcation(Yeff_7,R(6),R(15),L(6),L(15),RW1(15),RW2(15),CWK(15),E(6),E(15),h(6),h(15),v,v,rho,rho,omega(ih),a(6),a(15));

    %5-6-14 --> Type III Bifurcation
    [B5_A5,B14_A14,Yeff_5] = type_III_bifurcation(Yeff_6,R(5),R(14),L(5),L(14),RW1(14),RW2(14),CWK(14),E(5),E(14),h(5),h(14),v,v,rho,rho,omega(ih),a(5),a(14));
    
    %4-5 --> Special case of two connected vessels
    s4out = (R(4)-tan(a(4))*L(4))/sin(a(4));
    [J_13s4out,Y_13s4out,J_43s4out,Y_43s4out,fs4out] = besselfunctions(a(4),s4out,omega(ih),E(4),rho,v,h(4));      
    Y_s4out = (2*pi*(1-cos(a(4))))*(fs4out/rho)^0.5*s4out^2.5;
    B4_A4 = -(1i*Y_s4out*J_43s4out+Yeff_5*J_13s4out)/(Yeff_5*Y_13s4out+1i*Y_s4out*Y_43s4out);
    
    s4in = (R(4)-tan(a(4))*0)/sin(a(4));
    [J_13s4in,Y_13s4in,J_43s4in,Y_43s4in,fs4in] = besselfunctions(a(4),s4in,omega(ih),E(4),rho,v,h(4));      
    Y_s4in = (2*pi*(1-cos(a(4))))*(fs4in/rho)^0.5*s4in^2.5;
    Yeff_4 = -1i*Y_s4in*(J_43s4in+B4_A4*Y_43s4in)/(J_13s4in+B4_A4*Y_13s4in);
    
    %3-4-13 --> Type III Bifurcation  
    [B3_A3,B13_A13,Yeff_3] = type_III_bifurcation(Yeff_4,R(3),R(13),L(3),L(13),RW1(13),RW2(13),CWK(13),E(3),E(13),h(3),h(13),v,v,rho,rho,omega(ih),a(3),a(13));
 
    %2-3-12 --> Type III Bifurcation    
    [B2_A2,B12_A12,Yeff_2] = type_III_bifurcation(Yeff_3,R(2),R(12),L(2),L(12),RW1(12),RW2(12),CWK(12),E(2),E(12),h(2),h(12),v,v,rho,rho,omega(ih),a(2),a(12));
    
    %1-2-11 --> Type V Bifurcation    
    [A1,B1,B11_A11] = type_V_bifurcation(Yeff_2,R(1),R(11),L(1),L(11),RW1(11),RW2(11),CWK(11),E(1),E(11),h(1),h(11),v,v,rho,rho,omega(ih),a(1),a(11));
    
    %% Forward Propagation 
    %Inlet of vessel 1
    x1 = 0;
    s1 = (R(1)-tan(a(1))*x1)/sin(a(1));
    [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a(1),s1,omega(ih),E(1),rho,v,h(1)); 
    Y_s1 = (2*pi*(1-cos(a(1))))*(fs1/rho)^0.5*s1^2.5;
    Q1(ih) = -(1i*Y_s1*(s1^-0.5)*(A1*J_43s1+B1*Y_43s1));
    P1(ih) = ((s1^-0.5)*(A1*J_13s1+B1*Y_13s1));
    
    %Outlet of vessel 1
    
    x1 = L(1);
    s1out = (R(1)-tan(a(1))*x1)/sin(a(1));
    [J_13s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a(1),s1out,omega(ih),E(1),rho,v,h(1)); 
    Y_s1out = (2*pi*(1-cos(a(1))))*((fs1out/rho)^0.5)*s1out^2.5;
    Q1out(ih) = -(1i*Y_s1out*(s1out^-0.5)*(A1*J_43s1out+B1*Y_43s1out));
    P1out(ih) = ((s1out^-0.5)*(A1*J_13s1out+B1*Y_13s1out));
    
    % Brachiocephalic

    [Q11out(ih), P11out(ih),A,B] = vessel(P1out(ih),L(11),R(11),a(11),omega(ih),E(11),rho,v,h(11),B11_A11);
    
    %AO II
    
    [Q2out(ih), P2out(ih)] = vessel(P1out(ih),L(2),R(2),a(2),omega(ih),E(2),rho,v,h(2),B2_A2);

    %L com. carotid
    
    [Q12out(ih), P12out(ih)] = vessel(P2out(ih),L(12),R(12),a(12),omega(ih),E(12),rho,v,h(12),B12_A12);

    %AO III

    [Q3out(ih), P3out(ih)] = vessel(P2out(ih),L(3),R(3),a(3),omega(ih),E(3),rho,v,h(3),B3_A3);

    %Left subclavian

    [Q13out(ih), P13out(ih)] = vessel(P3out(ih),L(13),R(13),a(13),omega(ih),E(13),rho,v,h(13),B13_A13);

    %AO IV

    [Q4out(ih), P4out(ih)] = vessel(P3out(ih),L(4),R(4),a(4),omega(ih),E(4),rho,v,h(4),B4_A4);

    %AO V   

    [Q5out(ih), P5out(ih)] = vessel(P4out(ih),L(5),R(5),a(5),omega(ih),E(5),rho,v,h(5),B5_A5);

    %AO VI     

    [Q6out(ih), P6out(ih)] = vessel(P5out(ih),L(6),R(6),a(6),omega(ih),E(6),rho,v,h(6),B6_A6);

    %Sup. mesenteric

    [Q15out(ih), P15out(ih)] = vessel(P6out(ih),L(15),R(15),a(15),omega(ih),E(15),rho,v,h(15),B15_A15);

    %AO VII
   
    [Q7out(ih), P7out(ih)] = vessel(P6out(ih),L(7),R(7),a(7),omega(ih),E(7),rho,v,h(7),B7_A7);

    %Renal
    
    [Q16out(ih), P16out(ih)] = vessel(P7out(ih),L(16),R(16),a(16),omega(ih),E(16),rho,v,h(16),B16_A16);

    %AO VIII
    
    [Q8out(ih), P8out(ih)] = vessel(P7out(ih),L(8),R(8),a(8),omega(ih),E(8),rho,v,h(8),B8_A8);

    %AO IX

    [Q9out(ih), P9out(ih)] = vessel(P8out(ih),L(9),R(9),a(9),omega(ih),E(9),rho,v,h(9),B9_A9);

    %Inf. mesenteric

    [Q18out(ih), P18out(ih)] = vessel(P9out(ih),L(18),R(18),a(18),omega(ih),E(18),rho,v,h(18),B18_A18);

    %AO X

    [Q10out(ih), P10out(ih)] = vessel(P9out(ih),L(10),R(10),a(10),omega(ih),E(10),rho,v,h(10),B10_A10);

    %R com. iliac
    [Q19out(ih), P19out(ih)] = vessel(P10out(ih),L(19),R(19),a(19),omega(ih),E(19),rho,v,h(19),B19_A19);

end

%% Inverse Fourier
q1 = zeros(size(t)) + real(Q1(1)*F(1));
p1 = zeros(size(t)) + real(P1(1)*F(1));
q11 = zeros(size(t)) + real(Q11out(1)*F(1));
p11 = zeros(size(t)) + real(P11out(1)*F(1));
q12 = zeros(size(t)) + real(Q12out(1)*F(1));
p12 = zeros(size(t)) + real(P12out(1)*F(1));
q13 = zeros(size(t)) + real(Q13out(1)*F(1));
p13 = zeros(size(t)) + real(P13out(1)*F(1));
q15 = zeros(size(t)) + real(Q15out(1)*F(1));
p15 = zeros(size(t)) + real(P15out(1)*F(1));
q16 = zeros(size(t)) + real(Q16out(1)*F(1));
p16 = zeros(size(t)) + real(P16out(1)*F(1));
q18 = zeros(size(t)) + real(Q18out(1)*F(1));
p18 = zeros(size(t)) + real(P18out(1)*F(1));
q19 = zeros(size(t)) + real(Q19out(1)*F(1));
p19 = zeros(size(t)) + real(P19out(1)*F(1));

for ih = 1:nh
     
        
      q1 = q1 + real(conj(Q1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p1 = p1 + real(conj(P1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p11 = p11 + real(conj(P11out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q11 = q11 + real(conj(Q11out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p12 = p12 + real(conj(P12out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q12 = q12 + real(conj(Q12out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p13 = p13 + real(conj(P13out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q13 = q13 + real(conj(Q13out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p15 = p15 + real(conj(P15out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q15 = q15 + real(conj(Q15out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p16 = p16 + real(conj(P16out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q16 = q16 + real(conj(Q16out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p18 = p18 + real(conj(P18out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q18 = q18 + real(conj(Q18out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p19 = p19 + real(conj(P19out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      q19 = q19 + real(conj(Q19out(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));


end

toc
%% Error calculations
%Average
yiliacoutletpressure = interp1q(riliac_outlet_19_Pressure(1:end,1),riliac_outlet_19_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    err19(i) = abs((yiliacoutletpressure(i) - p19(i))/(yiliacoutletpressure(i)))^2;

end
error19pressure = sqrt(mean(err19));

yiliacoutletflow = interp1q(riliac_outlet_19_flow(1:end,1),riliac_outlet_19_flow(1:end,2),xi);
for i = 2:length(t)-2
    
    err19(i) = abs((yiliacoutletflow(i)*10^-6 - q19(i))/(max(yiliacoutletflow)*10^-6))^2;

end
error19flow = sqrt(mean(err19));

yiimmaoutletpressure = interp1q(imma_outlet_18_Pressure(1:end,1),imma_outlet_18_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    err18(i) = abs((yiimmaoutletpressure(i) - p18(i))/(yiimmaoutletpressure(i)))^2;

end
error18pressure = sqrt(mean(err18));

yiimmaoutletflow = interp1q(imma_outlet_18_flow(1:end,1),imma_outlet_18_flow(1:end,2),xi);
for i = 2:length(t)-2
    
    err18(i) = abs((yiimmaoutletflow(i)*10^-6 - q18(i))/(max(yiimmaoutletflow)*10^-6))^2;

end
error18flow = sqrt(mean(err18));

yirenaoutletpressure = interp1q(rena_outlet_16_Pressure(1:end,1),rena_outlet_16_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    err16(i) = abs((yirenaoutletpressure(i) - p16(i))/(yirenaoutletpressure(i)))^2;

end
error16pressure = sqrt(mean(err16));

yirenaoutletflow = interp1q(rena_outlet_16_flow(1:end,1),rena_outlet_16_flow(1:end,2),xi);
for i = 2:length(t)-5
    
    err16(i) = ((yirenaoutletflow(i)*10^-6 - q16(i))/(max(yirenaoutletflow)*10^-6))^2;

end
error16flow = sqrt(mean(err16));

yismaoutletpressure = interp1q(sma_outlet_15_Pressure(1:end,1),sma_outlet_15_Pressure(1:end,2),xi);
for i = 2:length(t)-3
    
    err15(i) = abs((yismaoutletpressure(i) - p15(i))/(yismaoutletpressure(i)))^2;

end
error15pressure = sqrt(mean(err15));

yismaoutletflow = interp1q(sma_outlet_15_flow(1:end,1),sma_outlet_15_flow(1:end,2),xi);
for i = 2:length(t)-5
    
    err15(i) = abs((yismaoutletflow(i)*10^-6 - q15(i))/(max(yismaoutletflow)*10^-6))^2;

end
error15flow = sqrt(mean(err15));

yilsuboutletpressure = interp1q(lsub_outlet_13_Pressure(1:end,1),lsub_outlet_13_Pressure(1:end,2),xi);
for i = 2:length(t)-3
    
    err13(i) = abs((yilsuboutletpressure(i) - p13(i))/(yilsuboutletpressure(i)))^2;

end
error13pressure = sqrt(mean(err13));

yilsuboutletflow = interp1q(lsub_outlet_13_flow(1:end,1),lsub_outlet_13_flow(1:end,2),xi);
for i = 2:length(t)-4
    
    err13(i) = abs((yilsuboutletflow(i)*10^-6 - q13(i))/(max(yilsuboutletflow)*10^-6))^2;

end
error13flow = sqrt(abs(mean(err13)));

yilccaoutletpressure = interp1q(lcca_outlet_12_Pressure(1:end,1),lcca_outlet_12_Pressure(1:end,2),xi);
for i = 2:length(t)-3
    
    err12(i) = abs((yilccaoutletpressure(i) - p12(i))/(yilccaoutletpressure(i)))^2;

end
error12pressure = sqrt(mean(err12));

yilccaoutletflow = interp1q(lcca_outlet_12_flow(1:end,1),lcca_outlet_12_flow(1:end,2),xi);
for i = 2:length(t)-2
    
    err12(i) = ((yilccaoutletflow(i)*10^-6 - q12(i))/(max(yilccaoutletflow)*10^-6))^2;

end
error12flow = sqrt(abs(mean(err12)));

yibcoutletpressure = interp1q(bc_outlet_11_Pressure(1:end,1),bc_outlet_11_Pressure(1:end,2),xi);
for i = 2:length(t)-3
    
    err11(i) = abs((yibcoutletpressure(i) - p11(i))/(yibcoutletpressure(i)))^2;

end
error11pressure = sqrt(mean(err11));

yibcoutletflow = interp1q(bc_outlet_11_flow(1:end,1),bc_outlet_11_flow(1:end,2),xi);
for i = 3:length(t)-3
    
    err11(i) = ((yibcoutletflow(i)*10^-6 - q11(i))/(max(yibcoutletflow)*10^-6))^2;

end
error11flow = sqrt(mean(err11));

yiinletpressure = interp1q(Inlet_Pressure(1:end,1),Inlet_Pressure(1:end,2),xi);
for i = 2:length(t)-3
    
    err1(i) = abs((yiinletpressure(i) - p1(i))/(yiinletpressure(i)))^2;

end
error1pressure = sqrt(mean(err1));

%MAX
emax19flow = max(abs((q19(2:end-1)-yiliacoutletflow(2:end-1)*10^-6)/max(yiliacoutletflow(2:end-1)*10^-6)));
emax19pressure = max(abs((p19-yiliacoutletpressure)./yiliacoutletpressure));

emax18flow = max(abs((q18(2:end-1)-yiimmaoutletflow(2:end-1)*10^-6)/max(yiimmaoutletflow(2:end-1)*10^-6)));
emax18pressure = max(abs((p18-yiimmaoutletpressure)./yiimmaoutletpressure));

emax16flow = max(abs((q16(2:end-1)-yirenaoutletflow(2:end-1)*10^-6)/max(yirenaoutletflow(2:end-1)*10^-6)));
emax16pressure = max(abs((p16-yirenaoutletpressure)./yirenaoutletpressure));

emax15flow = max(abs((q15(2:end-1)-yismaoutletflow(2:end-1)*10^-6)/max(yismaoutletflow(2:end-1)*10^-6)));
emax15pressure = max(abs((p15-yismaoutletpressure)./yismaoutletpressure));

emax13flow = max(abs((q13(2:end-1)-yilsuboutletflow(2:end-1)*10^-6)/max(yilsuboutletflow(2:end-1)*10^-6)));
emax13pressure = max(abs((p13-yilsuboutletpressure)./yilsuboutletpressure));

emax12flow = max(abs((q12(2:end-1)-yilccaoutletflow(2:end-1)*10^-6)/max(yilccaoutletflow(2:end-1)*10^-6)));
emax12pressure = max(abs((p12-yilccaoutletpressure)./yilccaoutletpressure));

emax11flow = max(abs((q11(2:end-1)-yibcoutletflow(2:end-1)*10^-6)/max(yibcoutletflow(2:end-1)*10^-6)));
emax11pressure = max(abs((p11-yibcoutletpressure)./yibcoutletpressure));

emax1pressure = max(abs((p1-yiinletpressure)./yiinletpressure));

%SYS
esys19flow = (max(q19(2:end-1)-max(yiliacoutletflow(2:end-1))*10^-6))/max(yiliacoutletflow(2:end-1)*10^-6);
esys19pressure = (max(p19)-max(yiliacoutletpressure))/max(yiliacoutletpressure);

esys18flow = (max(q18(2:end-1))-max(yiimmaoutletflow(2:end-1))*10^-6)/(max(yiimmaoutletflow(2:end-1)*10^-6));
esys18pressure = (max(p18)-max(yiimmaoutletpressure))/max(yiimmaoutletpressure);

esys16flow = (max(q16(2:end-1))-max(yirenaoutletflow(2:end-1)*10^-6))/max(yirenaoutletflow(2:end-1)*10^-6);
esys16pressure = (max(p16)-max(yirenaoutletpressure))/max(yirenaoutletpressure);

esys15flow = (max(q15(2:end-1))-max(yismaoutletflow(2:end-1)*10^-6))/max(yismaoutletflow(2:end-1)*10^-6);
esys15pressure = (max(p15)-max(yismaoutletpressure))/max(yismaoutletpressure);

esys13flow = (max(q13(2:end-1))-max(yilsuboutletflow(2:end-1))*10^-6)/max(yilsuboutletflow(2:end-1)*10^-6);
esys13pressure = (max(p13)-max(yilsuboutletpressure))/max(yilsuboutletpressure);

esys12flow = (max(q12(2:end-1))-max(yilccaoutletflow(2:end-1))*10^-6)/max(yilccaoutletflow(2:end-1)*10^-6);
esys12pressure = (max(p12)-max(yilccaoutletpressure))/max(yilccaoutletpressure);

esys11flow = (max(q11(2:end-1))-max(yibcoutletflow(2:end-1)*10^-6))/max(yibcoutletflow(2:end-1)*10^-6);
esys11pressure = (max(p11)-max(yibcoutletpressure))/max(yibcoutletpressure);

esys1pressure = (max(p1)-max(yiinletpressure))/max(yiinletpressure);

%DIAS
edias19flow = (min(q19(2:end-1)-min(yiliacoutletflow(2:end-1))*10^-6))/max(yiliacoutletflow(2:end-1)*10^-6);
edias19pressure = (min(p19)-min(yiliacoutletpressure))/min(yiliacoutletpressure);

edias18flow = (min(q18(2:end-1))-min(yiimmaoutletflow(2:end-1))*10^-6)/(max(yiimmaoutletflow(2:end-1)*10^-6));
edias18pressure = (min(p18)-min(yiimmaoutletpressure))/min(yiimmaoutletpressure);

edias16flow = (min(q16(2:end-1))-min(yirenaoutletflow(2:end-1)*10^-6))/max(yirenaoutletflow(2:end-1)*10^-6);
edias16pressure = (min(p16)-min(yirenaoutletpressure))/min(yirenaoutletpressure);

edias15flow = (min(q15(2:end-1))-min(yismaoutletflow(2:end-1)*10^-6))/max(yismaoutletflow(2:end-1)*10^-6);
edias15pressure = (min(p15)-min(yismaoutletpressure))/min(yismaoutletpressure);

edias13flow = (min(q13(2:end-1))-min(yilsuboutletflow(2:end-1))*10^-6)/max(yilsuboutletflow(2:end-1)*10^-6);
edias13pressure = (min(p13)-min(yilsuboutletpressure))/min(yilsuboutletpressure);

edias12flow = (min(q12(2:end-1))-min(yilccaoutletflow(2:end-1))*10^-6)/max(yilccaoutletflow(2:end-1)*10^-6);
edias12pressure = (min(p12)-min(yilccaoutletpressure))/min(yilccaoutletpressure);

edias11flow = (min(q11(2:end-1))-min(yibcoutletflow(2:end-1)*10^-6))/max(yibcoutletflow(2:end-1)*10^-6);
edias11pressure = (min(p11)-min(yibcoutletpressure))/min(yibcoutletpressure);

edias1pressure = (min(p1)-min(yiinletpressure))/min(yiinletpressure);

errf = (error11flow+error12flow+error13flow+error15flow+error16flow+error18flow+error19flow)/7;
errp = (error11pressure+error12pressure+error13pressure+error15pressure+error16pressure+error18pressure+error19pressure+error1pressure)/8;

%% Plots
ld = 2.5;
colour1 = 'b';
colour2 = 'r';
fontxt = 21;
fontxt2 = 21;
figure
plot(t,q1*10^6,colour1,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = title('Volume flow rate - inlet BC (1)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;
axis([0 1.1,-100 600])

figure
plot(t,p1*10^-3,colour1,'LineWidth',ld)
hold on
plot(Inlet_Pressure(1:end,1),Inlet_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - inlet (1)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error1pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax1pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys1pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias1pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q11*10^6,colour1,'LineWidth',ld)
hold on
plot(bc_outlet_11_flow(1:end,1),bc_outlet_11_flow(1:end,2),colour2,'LineWidth',ld);
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - bc outlet (11)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = 13;
txt = ['avg%: ',num2str(round(error11flow*100,2))];
t1 = text(0.7,35,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax11flow*100,2))];
t1 = text(0.7,30,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys11flow*100,2))];
t1 = text(0.7,25,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias11flow*100,2))];
t1 = text(0.7,20,txt);
t1.FontSize = fontxt;
axis([0 1.1,-5 50])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p11*10^-3,colour1,'LineWidth',ld)
hold on
plot(bc_outlet_11_Pressure(1:end,1),bc_outlet_11_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure -  bc outlet (11)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error11pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax11pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys11pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias11pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q12*10^6,colour1,'LineWidth',ld)
hold on
plot(lcca_outlet_12_flow(1:end,1),lcca_outlet_12_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - lcca outlet (12)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error12flow*100,2))];
t1 = text(0.7,9.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax12flow*100,2))];
t1 = text(0.7,7.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys12flow*100,2))];
t1 = text(0.7,5.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias12flow*100,2))];
t1 = text(0.7,4.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,-5 15])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p12*10^-3,colour1,'LineWidth',ld)
hold on
plot(lcca_outlet_12_Pressure(1:end,1),lcca_outlet_12_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - lcca outlet (12)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error12pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax12pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys12pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias12pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q13*10^6,colour1,'LineWidth',ld)
hold on
plot(lsub_outlet_13_flow(1:end,1),lsub_outlet_13_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - lsub outlet (13)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error13flow*100,2))];
t1 = text(0.7,35,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax13flow*100,2))];
t1 = text(0.7,30,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys13flow*100,2))];
t1 = text(0.7,25,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias13flow*100,2))];
t1 = text(0.7,20,txt);
t1.FontSize = fontxt;
axis([0 1.1,-3 50])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p13*10^-3,colour1,'LineWidth',ld)
hold on
plot(lsub_outlet_13_Pressure(1:end,1),lsub_outlet_13_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - lsub outlet (13)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error13pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax13pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys13pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias13pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q15*10^6,colour1,'LineWidth',ld)
hold on
plot(sma_outlet_15_flow(1:end,1),sma_outlet_15_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - sma outlet (15)');t2.FontSize = 13;
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error15flow*100,2))];
t1 = text(0.7,35,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax15flow*100,2))];
t1 = text(0.7,30,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys15flow*100,2))];
t1 = text(0.7,25,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias15flow*100,2))];
t1 = text(0.7,20,txt);
t1.FontSize = fontxt;
axis([0 1.1,-3 50])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p15*10^-3,colour1,'LineWidth',ld)
hold on
plot(sma_outlet_15_Pressure(1:end,1),sma_outlet_15_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - sma outlet (15)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error15pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax15pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys15pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias15pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q16*10^6,colour1,'LineWidth',ld)
hold on
plot(rena_outlet_16_flow(1:end,1),rena_outlet_16_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - r ren (16)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error16flow*100,2))];
t1 = text(0.7,35,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax16flow*100,2))];
t1 = text(0.7,30,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys16flow*100,2))];
t1 = text(0.7,25,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias16flow*100,2))];
t1 = text(0.7,20,txt);
t1.FontSize = fontxt;
axis([0 1.1,-3 50])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p16*10^-3,colour1,'LineWidth',ld)
hold on
plot(rena_outlet_16_Pressure(1:end,1),rena_outlet_16_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - r ren (16)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error16pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax16pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys16pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias16pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q18*10^6,colour1,'LineWidth',ld)
hold on
plot(imma_outlet_18_flow(1:end,1),imma_outlet_18_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - imma outlet (18)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error18flow*100,2))];
t1 = text(0.7,9.2,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax18flow*100,2))];
t1 = text(0.7,7.5,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys18flow*100,2))];
t1 = text(0.7,5.8,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias18flow*100,2))];
t1 = text(0.7,4.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,-5 15])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p18*10^-3,colour1,'LineWidth',ld)
hold on
plot(imma_outlet_18_Pressure(1:end,1),imma_outlet_18_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - imma outlet (18)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error18pressure*100,2))];
t1 = text(0.7,14.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax18pressure*100,2))];
t1 = text(0.7,13.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys18pressure*100,2))];
t1 = text(0.7,12.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias18pressure*100,2))];
t1 = text(0.7,12.1,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 17])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,q19*10^6,colour1,'LineWidth',ld)
hold on
plot(riliac_outlet_19_flow(1:end,1),riliac_outlet_19_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - r iliac outlet (19)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error19flow*100,2))];
t1 = text(0.7,45,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax19flow*100,2))];
t1 = text(0.7,37,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys19flow*100,2))];
t1 = text(0.7,29,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias19flow*100,2))];
t1 = text(0.7,21,txt);
t1.FontSize = fontxt;
axis([0 1.1,-25 70])
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

figure
plot(t,p19*10^-3,colour1,'LineWidth',ld)
hold on
plot(riliac_outlet_19_Pressure(1:end,1),riliac_outlet_19_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - r iliac outlet (19)');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(error19pressure*100,2))];
t1 = text(0.7,16,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emax19pressure*100,2))];
t1 = text(0.7,15,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esys19pressure*100,2))];
t1 = text(0.7,14,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(edias19pressure*100,2))];
t1 = text(0.7,13,txt);
t1.FontSize = fontxt;
axis([0 1.1,8 19]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

