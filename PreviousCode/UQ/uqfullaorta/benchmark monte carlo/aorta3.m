%% Function for implementation of the full aorta model
function [Qout7,l,t] = aorta3(R1)
%% Flores et all (2016), Table 5, Parameters of the full-aorta model
L = [7.0357,0.8,0.9,6.4737,15.2,1.8,0.7,0.7,4.3,4.3,3.4,3.4,3.4,3.2,6,3.2,3.2,5,8.5,8.5]*10^-2;
Rin = R1;
h = 0.1*[15.2,13.9,13.7,13.5,12.3,9.9,9.7,9.62,9.55,9.07,6.35,3.6,4.8,4.45,3.75,2.8,2.8,2,6,6]*10^-3; %h was chosen to be 10% of Rd
%Rin = Rin;
Rout = [Rin(2)*10^3,Rin(3)*10^3,Rin(4)*10^3,Rin(5)*10^3,Rin(6)*10^3,Rin(7)*10^3,Rin(8)*10^3,Rin(9)*10^3,Rin(10)*10^3,8.6,Rin(11)*10^3,Rin(12)*10^3,Rin(13)*10^3,Rin(14)*10^3,Rin(15)*10^3,Rin(16)*10^3,Rin(17)*10^3,Rin(18)*10^3,Rin(19)*10^3,Rin(20)*10^3]*10^-3;

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
load('/Users/georgiosefstathiou/Desktop/MSc/Project/Code/automeris/aorta/BC.csv')

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
%nh = N/2;
nh = 10;
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
for ih = 0:nh   %Work with the first 10 harmonics as in Papadakis 2019     
    
    %% Freequency
    
    if ih == 0
         omega(ih+1) = 10^-10;
    
    else
         omega(ih+1) = -ih*om;
    end
    
    ih=ih+1;
   
    %% Backward Propagation  

    %10-19-20 --> Type II Bifurcation
    [B10_A10,B19_A19,B20_A20,Yeff_10] = type_II_bifurcation(R(10),R(19),R(20),L(10),L(19),L(20),RW1(19),RW2(19),CWK(19),RW1(20),RW2(20),CWK(20),E(10),E(19),E(20),h(10),h(19),h(20),v,v,v,rho,rho,rho,omega(ih),a(10),a(19),a(20));
%     
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
%p1 = zeros(size(t)) + real(P1(1)*F(1));
q11 = zeros(size(t)) + real(Q11out(1)*F(1));
%p11 = zeros(size(t)) + real(P11out(1)*F(1));
q12 = zeros(size(t)) + real(Q12out(1)*F(1));
%p12 = zeros(size(t)) + real(P12out(1)*F(1));
q13 = zeros(size(t)) + real(Q13out(1)*F(1));
%p13 = zeros(size(t)) + real(P13out(1)*F(1));
q15 = zeros(size(t)) + real(Q15out(1)*F(1));
%p15 = zeros(size(t)) + real(P15out(1)*F(1));
q16 = zeros(size(t)) + real(Q16out(1)*F(1));
%p16 = zeros(size(t)) + real(P16out(1)*F(1));
q18 = zeros(size(t)) + real(Q18out(1)*F(1));
%p18 = zeros(size(t)) + real(P18out(1)*F(1));
q19 = zeros(size(t)) + real(Q19out(1)*F(1));
%p19 = zeros(size(t)) + real(P19out(1)*F(1));

for ih = 1:nh
     
        
      q1 = q1 + real(Q1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
    %  p1 = p1 + real(P1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
    %  p11 = p11 + real(P11out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      q11 = q11 + real(Q11out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
     % p12 = p12 + real(P12out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      q12 = q12 + real(Q12out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
     % p13 = p13 + real(P13out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      q13 = q13 + real(Q13out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
     % p15 = p15 + real(P15out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      q15 = q15 + real(Q15out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
     % p16 = p16 + real(P16out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      q16 = q16 + real(Q16out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
     % p18 = p18 + real(P18out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      q18 = q18 + real(Q18out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
    %  p19 = p19 + real(P19out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
    q19 = q19 + real(Q19out(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));


end

l = length(q1);
% Qout1 = q11;
% Qout2 = q12;
% Qout3 = q13;
% Qout4 = q15;
% Qout5 = q16;
% Qout6 = q18;
Qout7 = q19;
