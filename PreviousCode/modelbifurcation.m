%% Programme for implementation of the model bifurcation case
clear;
clc;

%% Defining the constant parameters of the problem

%Angle a
a = 0.001;
a2 = 0.001;

%Parent vessel (1)
E1 = 500*10^3;                                       %Table 3 (Flores 2016)
rho1 = 1060;                                         %Table 3 (Flores 2016)
v1 = 0.5;                                            %Table 3 (Flores 2016)
h1 = 1.032*10^-3;                                    %Table 3 (Flores 2016)
R1 = 0.89*10^-2;
L1 = 8.6*10^-2;
c1 = ((E1*h1)/((1-v1^2)*rho1*2*R1))^0.5;
Y1 = (pi*R1^2)/(rho1*c1);

%Daughter vessels (2) and (3)
E2 = 700*10^3;                                       %Table 3 (Flores 2016)
rho2 = 1060;                                         %Table 3 (Flores 2016)
v2 = 0.5;                                            %Table 3 (Flores 2016)
h2 = 0.72*10^-3;                                     %Table 3 (Flores 2016)
R2 = 0.6125*10^-2;
L2 = 8.5*10^-2;
c2 = ((E2*h2)/((1-v2^2)*rho2*2*R2))^0.5;
Y2 = (pi*R2^2)/(rho2*c2);
Y3 = Y2;


E3 = 700*10^3;                                       %Table 3 (Flores 2016)
rho3 = 1060;                                         %Table 3 (Flores 2016)
v3 = 0.5;                                            %Table 3 (Flores 2016)
h3 = 0.72*10^-3;                                     %Table 3 (Flores 2016)
R3 = 0.6125*10^-2;
L3 = 8.5*10^-2;

%Impedence calculation data
RW1 = 6.8123*10^7;
RW2 = 3.1013*10^9;
Cwk = 3.6664*10^-10;

RW13 = 6.8123*10^7;
RW23 = 3.1013*10^9;
Cwk3 = 3.6664*10^-10;

% Optimised Values

% E1 = 500*10^3;                                       %Table 3 (Flores 2016)
% rho1 = 1060;                                         %Table 3 (Flores 2016)
% v1 = 0.5;                                            %Table 3 (Flores 2016)
% be1 = 0.00107412679934682;
% h1 = 0.75/(E1*be1);                                    %Table 3 (Flores 2016)
% R1 = 0.00895572408393668;
% L1 = 0.105806045268162;
% c1 = ((E1*h1)/((1-v1^2)*rho1*2*R1))^0.5;
% Y1 = (pi*R1^2)/(rho1*c1);
% 
% %Daughter vessels (2) and (3)
% E2 = 700*10^3;                                       %Table 3 (Flores 2016)
% rho2 = 1060;                                         %Table 3 (Flores 2016)
% v2 = 0.5;                                            %Table 3 (Flores 2016)
% be2 = 0.00104673896585901;
% h2 = 0.75/(E2*be2);                                     %Table 3 (Flores 2016)
% R2 = 0.00592762150465316;
% L2 = 0.133512189773607;
% c2 = ((E2*h2)/((1-v2^2)*rho2*2*R2))^0.5;
% Y2 = (pi*R2^2)/(rho2*c2);
% Y3 = Y2;
% 
% 
% E3 = E2;                                       %Table 3 (Flores 2016)
% rho3 = 1060;                                         %Table 3 (Flores 2016)
% v3 = 0.5;                                            %Table 3 (Flores 2016)
% h3 = h2;                                     %Table 3 (Flores 2016)
% R3 = R2;
% L3 = L2;
% 
% %Impedence calculation data
% RW1 = 86210920.7522733;
% RW2 = 3077962905.34146;
% Cwk = 3.54153142218048e-10;
% 
% RW13 = RW1;
% RW23 = RW2;
% Cwk3 = Cwk;

%% Importing the data from the Flores plots
load('grabit/ao_mid_flow.mat')
load('automeris/bifurcation/Inlet_ao_Pressure.csv')
load('automeris/bifurcation/mid_ao_Pressure.csv')
load('automeris/bifurcation/junction_Pressure.csv')
load('automeris/bifurcation/mid_il_Pressure.csv')
load('automeris/bifurcation/outlet_il_Pressure.csv')
load('automeris/bifurcation/BC.csv')
load('automeris/bifurcation/junction_flow.csv')
load('automeris/bifurcation/mid_il_flow.csv')
load('automeris/bifurcation/outlet_il_flow.csv')


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

%% Frequency domain calculations
for ih = 0:nh   
%% Backward Propagation       
     
    if ih == 0
        omega(ih+1) = 10^-10;
    
    else
         omega(ih+1) = -ih*om;
    end
    
    ih=ih+1;
    
    %% Step 1
    %Firstly the terms B2/A2 and B3/A3 need to be calculated
    
    %Impedance Z (same for both daughter vessels)
    Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);
    Z3 = (RW13+RW23-1i*omega(ih)*RW13*RW23*Cwk3)/(1-1i*omega(ih)*RW23*Cwk3);

    %Bessel functions at s2out and s3out
    s2out = (R2-tan(a2)*L2)/sin(a2);
    s3out = (R3-tan(a2)*L3)/sin(a2);
    [J_13_s2out,Y_13s2out,J_43s2out,Y_43s2out,fs2out] = besselfunctions(a2,s2out,omega(ih),E2,rho2,v2,h2); 
    [J_13_s3out,Y_13s3out,J_43s3out,Y_43s3out,fs3out] = besselfunctions(a2,s3out,omega(ih),E3,rho3,v3,h3); 
    
    Y_s2out = (2*pi*(1-cos(a2)))*(fs2out/rho2)^0.5*s2out^2.5;
    Y_s3out = (2*pi*(1-cos(a2)))*(fs3out/rho2)^0.5*s3out^2.5;
    B2_A2 = -(J_13_s2out+1i*Y_s2out*J_43s2out*Z)/(Y_13s2out+1i*Y_s2out*Y_43s2out*Z);
    B3_A3 = -(J_13_s3out+1i*Y_s3out*J_43s3out*Z3)/(Y_13s3out+1i*Y_s3out*Y_43s3out*Z3);
    
    
    %% Step 2
    %Calculation of Yeff(s2in) and Yeff(s3in)
    s2in = R2/sin(a2);
    s3in = R3/sin(a2);

    [J_13_s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),E2,rho2,v2,h2); 
    [J_13_s3in,Y_13s3in,J_43s3in,Y_43s3in,fs3in] = besselfunctions(a2,s3in,omega(ih),E3,rho3,v3,h3); 

    Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
    Y_s3in = (2*pi*(1-cos(a2)))*(fs3in/rho3)^0.5*s3in^2.5;

    Yeff_s2in = -1i*Y_s2in*(J_43s2in+B2_A2*Y_43s2in)/(J_13_s2in+B2_A2*Y_13s2in);
    Yeff_s3in = -1i*Y_s3in*(J_43s3in+B3_A3*Y_43s3in)/(J_13_s3in+B3_A3*Y_13s3in);

    
    %% Step 3
    %Calculation of B1/A1
    s1out = (R1-tan(a)*L1)/sin(a);
    [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),E1,rho1,v1,h1);
    Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho1)^0.5*s1out^2.5;
    B1_A1  = -((Yeff_s2in+Yeff_s3in)*J_13_s1out+1i*Y_s1out*J_43s1out)/(1i*Y_s1out*Y_43s1out+(Yeff_s2in+Yeff_s3in)*Y_13s1out);
    
    %% Calculation of flow rate at a specific point
    x0 = 0;    %Inlet
    s0 = (R1-tan(a)*x0)/sin(a);
    [J_13_s0,Y_13s0,J_43s0,Y_43s0,fs0] = besselfunctions(a,s0,omega(ih),E1,rho1,v1,h1);
    Y_s0 = (2*pi*(1-cos(a)))*(fs0/rho1)^0.5*s0^2.5;
    
    x1 = 0;    %Inlet of parent vessel
    s1 = (R1-tan(a)*x1)/sin(a);
    [J_13_s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E1,rho1,v1,h1);
    Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho1)^0.5*s1^2.5; 
    Q1(ih) = (Y_s1*(s1^-0.5)*(J_43s1+B1_A1*Y_43s1))/(Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
    P1(ih) = -((s1^-0.5)*((J_13_s1)+B1_A1*Y_13s1))/(1i*Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
    
    x1 = L1/2;    %Middle of parent vessel
    s1 = (R1-tan(a)*x1)/sin(a);
    [J_13_s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E1,rho1,v1,h1);
    Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho1)^0.5*s1^2.5;
    Q1mid(ih) = (Y_s1*(s1^-0.5)*(J_43s1+B1_A1*Y_43s1))/(Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
    P1mid(ih) = -((s1^-0.5)*((J_13_s1)+B1_A1*Y_13s1))/(1i*Y_s0*(s0^-0.5)*(J_43s0+B1_A1*Y_43s0));
    
    x1 = L1;    %Junction
    s1out = (R1-tan(a)*x1)/sin(a);
    [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),E1,rho1,v1,h1);
    Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho1)^0.5*s1out^2.5;
    A1tilda = -1/(1i*Y_s0*(s0^(-0.5)*(J_43s0+B1_A1*Y_43s0)));
    B1tilda = B1_A1*A1tilda;
    P1junc(ih) = A1tilda*(s1out^-0.5)*((J_13_s1out)+B1_A1*Y_13s1out);
    Q1junc(ih) = -1i*(Y_s1out*(s1out^-0.5))*(A1tilda*J_43s1out+B1tilda*Y_43s1out);

%% Forward Propagation
    
   %% Step 1
    s2in = R2/sin(a2);
    [J_13s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),E2,rho2,v2,h2);
    Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
    A2tilda =  Q1junc(ih)/(2*(-1i*Y_s2in*(s2in^(-0.5))*(J_43s2in+B2_A2*Y_43s2in)));
    B2tilda = B2_A2*A2tilda;
   
   %% Step 2
   x2in = 0;    %Inlet at the daughter vessel
   s2in = (R2-tan(a2)*x2in)/sin(a2);
   [J_13s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),E2,rho2,v2,h2);
   Y_s2in = (2*pi*(1-cos(a2)))*((fs2in/rho2)^0.5)*s2in^2.5;
   Q2in(ih) = A2tilda*(-1i*Y_s2in*(s2in^(-0.5))*(J_43s2in+B2_A2*Y_43s2in));
   P2in(ih) = A2tilda*((s2in^(-0.5))*(J_13s2in+B2_A2*Y_13s2in));
   
   x2 = L2/2;    %Midlle of the daughter vessel
   s2 = (R2-tan(a2)*x2)/sin(a2);
   [J_13s2,Y_13s2,J_43s2,Y_43s2,fs2] = besselfunctions(a2,s2,omega(ih),E2,rho2,v2,h2);
   Y_s2 = (2*pi*(1-cos(a2)))*((fs2/rho2)^0.5)*s2^2.5;
   Q2ilmid(ih) = A2tilda*(-1i*Y_s2*(s2^(-0.5))*(J_43s2+B2_A2*Y_43s2));  
   P2ilmid(ih) = A2tilda*((s2^(-0.5))*(J_13s2+B2_A2*Y_13s2));
   
   x2 = L2;     %Outlet of the daughter vessel
   s2 = (R2-tan(a2)*x2)/sin(a2);
   [J_13_s2,Y_13s2,J_43s2,Y_43s2,fs2] = besselfunctions(a2,s2,omega(ih),E2,rho2,v2,h2);
   Y_s2 = (2*pi*(1-cos(a2)))*((fs2/rho2)^0.5)*s2^2.5;
   Q2ilout(ih) = A2tilda*(-1i*Y_s2*(s2^(-0.5))*(J_43s2+B2_A2*Y_43s2));
   P2ilout(ih) = Q2ilout(ih)*Z;
   
end

%% Inverse Fourier
q1 = zeros(size(t)) + real(F(1));
q1mid = zeros(size(t)) + real(Q1mid(1)*F(1));
q1junc = zeros(size(t)) + real(Q1junc(1)*F(1));
q2ilmid = zeros(size(t)) + real(Q2ilmid(1)*F(1));
q2ilout = zeros(size(t)) + real(Q2ilout(1)*F(1));
p1in = zeros(size(t)) + real(P1(1)*F(1));
p1mid = zeros(size(t)) + real(P1mid(1)*F(1));
p1junc = zeros(size(t)) + real(P1junc(1)*F(1));
p2in = zeros(size(t)) + real(P2in(1)*F(1));
p2mid = zeros(size(t)) + real(P2ilmid(1)*F(1));
p2ilout = zeros(size(t)) + real(P2ilout(1)*F(1));


 for ih = 1:nh
     
       q1 = q1 + real(Q1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       q1mid = q1mid + real(Q1mid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       q1junc = q1junc + real(Q1junc(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       q2ilmid = q2ilmid + real(Q2ilmid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       q2ilout = q2ilout + real(Q2ilout(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       p1in = p1in + real(P1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       p1mid = p1mid + real(P1mid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       p1junc = p1junc + real(P1junc(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       p2in = p2in + real(P2in(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       p2mid = p2mid + real(P2ilmid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
       p2ilout = p2ilout + real(P2ilout(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
   
 end




%% Error calculations
%Average
yiinletpressure = interp1q(Inlet_ao_Pressure(1:end,1),Inlet_ao_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errinletp(i) = abs((yiinletpressure(i) - p1in(i))/(yiinletpressure(i)))^2;

end
errorinletpressure = sqrt(mean(errinletp));


yiaomidpressure = interp1q(mid_ao_Pressure(1:end,1),mid_ao_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errmidp(i) = abs((yiaomidpressure(i) - p1mid(i))/(yiaomidpressure(i)))^2;

end
erroraomidpressure = sqrt(mean(errmidp));

yiaomidflow = interp1q(ao_mid_flow(1:end,1),ao_mid_flow(1:end,2),xi);
for i = 2:length(t)-2
    
     errmidf(i) = abs((yiaomidflow(i)*10^-6 - q1mid(i))/(max(yiaomidflow)*10^-6))^2;

end
erroraomidflow = sqrt(mean(errmidf));

yijuncpressure = interp1q(junction_Pressure(1:end,1),junction_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errjuncp(i) = abs((yijuncpressure(i) - p1junc(i))/(yijuncpressure(i)))^2;

end
errorjuncpressure = sqrt(mean(errjuncp));

yijuncflow = interp1q(junction_flow(1:end,1),junction_flow(1:end,2),xi);
for i = 1:length(t)-12
    
     errjuncf(i) = abs((yijuncflow(i)*10^-6 - q1junc(i))/(max(yijuncflow)*10^-6))^2;

end
errorjuncflow = sqrt(mean(errjuncf));

yiilmidpressure = interp1q(mid_il_Pressure(1:end,1),mid_il_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errillmidp(i) = abs((yiilmidpressure(i) - p2mid(i))/(yiilmidpressure(i)))^2;

end
errorilmidpressure = sqrt(mean(errillmidp));

yiilmidflow = interp1q(mid_il_flow(1:end,1),mid_il_flow(1:end,2),xi);
for i = 3:length(t)-3
    
     errilmidf(i) = abs((yiilmidflow(i)*10^-6 - q2ilmid(i))/(max(yiilmidflow)*10^-6))^2;

end
errorilmidflow = sqrt(mean(errilmidf));

yiiloutpressure = interp1q(outlet_il_Pressure(1:end,1),outlet_il_Pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    erriloutp(i) = abs((yiiloutpressure(i) - p2ilout(i))/(yiiloutpressure(i)))^2;

end
erroriloutpressure = sqrt(mean(erriloutp));

yiiloutflow = interp1q(outlet_il_flow(1:end,1),outlet_il_flow(1:end,2),xi);
for i = 7:length(t)-10
    
     erriloutf(i) = abs((yiiloutflow(i)*10^-6 - q2ilout(i))/(max(yiiloutflow)*10^-6))^2;

end
erroriloutflow = sqrt(mean(erriloutf));

%MAX
emaxaomidflow = max(abs((q1mid(2:end-1)-yiaomidflow(2:end-1)*10^-6)/max(yiaomidflow(2:end-1)*10^-6)));
emaxilmidflow = max(abs((q2ilmid(2:end-1)-yiilmidflow(2:end-1)*10^-6)/max(yiilmidflow(2:end-1)*10^-6)));
emaxjuncflow = max(abs((q1junc(2:end-1)-yijuncflow(2:end-1)*10^-6)/max(yijuncflow(2:end-1)*10^-6)));
emaxoutflow = max(abs((q2ilout(2:end-1)-yiiloutflow(2:end-1)*10^-6)/max(yiiloutflow(2:end-1)*10^-6)));
emaxinletpressure = max(abs((p1in(2:end-1)-yiinletpressure(2:end-1))/max(yiinletpressure(2:end-1))));
emaxaomidpressure = max(abs((p1mid(2:end-1)-yiaomidpressure(2:end-1))/max(yiaomidpressure(2:end-1))));
emaxjuncpressure = max(abs((p1junc(2:end-1)-yijuncpressure(2:end-1))/max(yijuncpressure(2:end-1))));
emaxilmidpressure = max(abs((p2mid(2:end-1)-yiilmidpressure(2:end-1))/max(yiilmidpressure(2:end-1))));
emaxoutpressure = max(abs((p2ilout(2:end-1)-yiiloutpressure(2:end-1))/max(yiiloutpressure(2:end-1))));


%SYS
esysaomidflow = (max(q1mid(2:end-1)-max(yiaomidflow(2:end-1))*10^-6))/max(yiaomidflow(2:end-1)*10^-6);
esysilmidflow = (max(q2ilmid(2:end-1)-max(yiilmidflow(2:end-1))*10^-6))/max(yiilmidflow(2:end-1)*10^-6);
esysjuncflow = (max(q1junc(2:end-1)-max(yijuncflow(2:end-1))*10^-6))/max(yijuncflow(2:end-1)*10^-6);
esysoutflow = (max(q2ilout(2:end-1)-max(yiiloutflow(2:end-1))*10^-6))/max(yiiloutflow(2:end-1)*10^-6);
esysinletpressure = (max(p1in(2:end-1)-max(yiinletpressure(2:end-1))))/max(yiinletpressure(2:end-1));
esysaomidpressure = (max(p1mid(2:end-1)-max(yiaomidpressure(2:end-1))))/max(yiaomidpressure(2:end-1));
esysjuncpressure = (max(p1junc(2:end-1)-max(yijuncpressure(2:end-1))))/max(yijuncpressure(2:end-1));
esysilmidpressure = (max(p2mid(2:end-1)-max(yiilmidpressure(2:end-1))))/max(yiilmidpressure(2:end-1));
esysoutpressure = (max(p2ilout(2:end-1)-max(yiiloutpressure(2:end-1))))/max(yiiloutpressure(2:end-1));

%DIAS
ediasaomidflow = (min(q1mid(2:end-1)-min(yiaomidflow(2:end-1))*10^-6))/max(yiaomidflow(2:end-1)*10^-6);
ediasilmidflow = (min(q2ilmid(2:end-1)-min(yiilmidflow(2:end-1))*10^-6))/max(yiilmidflow(2:end-1)*10^-6);
ediasjuncflow = (min(q1junc(2:end-1)-min(yijuncflow(2:end-1))*10^-6))/max(yijuncflow(2:end-1)*10^-6);
ediasoutflow = (min(q2ilout(2:end-1)-min(yiiloutflow(2:end-1))*10^-6))/max(yiiloutflow(2:end-1)*10^-6);
ediasinletpressure = (min(p1in(2:end-1)-min(yiinletpressure(2:end-1))))/max(yiinletpressure(2:end-1));
ediasaomidpressure = (min(p1mid(2:end-1)-min(yiaomidpressure(2:end-1))))/max(yiaomidpressure(2:end-1));
ediasjuncpressure = (min(p1junc(2:end-1)-min(yijuncpressure(2:end-1))))/max(yijuncpressure(2:end-1));
ediasilmidpressure = (min(p2mid(2:end-1)-min(yiilmidpressure(2:end-1))))/max(yiilmidpressure(2:end-1));
ediasoutpressure = (min(p2ilout(2:end-1)-min(yiiloutpressure(2:end-1))))/max(yiiloutpressure(2:end-1));

%% Ploting the results
ld = 2.5;
colour1 = 'b';
colour2 = 'r';
fontxt = 21;
fontxt2 = 21;

%Flow
figure
plot(t,Qin*10^6,colour1,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = title('Volume flow rate - inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time(s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q(ml/s)');
t2.FontSize = fontxt2;
axis([0 1.2,-40 100]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

% figure
% plot(t,q1mid*10^6,colour1,'LineWidth',ld)
% hold on
% plot(ao_mid_flow(1:end,1),ao_mid_flow(1:end,2),colour2,'LineWidth',ld)
% grid on
% t2 = legend('in-house','3-D');
% t2.FontSize = fontxt2;
% t2 = title('Volume flow rate - middle of the aorta');
% t2.FontSize = fontxt2;
% t2 = xlabel('time(s)');
% t2.FontSize = fontxt2;
% t2 = ylabel('Q(ml/s)');
% t2.FontSize = fontxt2;
% txt = ['avg%: ',num2str(round(erroraomidflow*100,2))];
% t1 = text(0.7,60,txt);
% t1.FontSize = fontxt;
% txt = ['max%: ',num2str(round(emaxaomidflow*100,2))];
% t1 = text(0.7,47,txt);
% t1.FontSize = fontxt;
% txt = ['sys%: ',num2str(round(esysaomidflow*100,2))];
% t1 = text(0.7,34,txt);
% t1.FontSize = fontxt;
% txt = ['dias%: ',num2str(round(ediasaomidflow*100,2))];
% t1 = text(0.7,21,txt);
% t1.FontSize = fontxt;
% axis([0 1.2,-40 100]);
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 0.5;
% ax.Layer = 'top';
% ax.XAxis.FontSize = fontxt2;
% ax.YAxis.FontSize = fontxt2;

figure
plot(t,q1junc*10^6,colour1,'LineWidth',ld)
hold on
plot(junction_flow(1:end,1),junction_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
txt = ['avg%: ',num2str(round(errorjuncflow*100,2))];
t1 = text(0.7,60,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxjuncflow*100,2))];
t1 = text(0.7,47,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysjuncflow*100,2))];
t1 = text(0.7,34,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasjuncflow*100,2))];
t1 = text(0.7,21,txt);
t1.FontSize = fontxt;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - junction');
t2.FontSize = fontxt2;
t2 = xlabel('time(s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q(ml/s)');
t2.FontSize = fontxt2;
axis([0 1.2,-40 100]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

% figure
% plot(t,q2ilmid*10^6,colour1,'LineWidth',ld)
% hold on
% plot(mid_il_flow(1:end,1),mid_il_flow(1:end,2),colour2,'LineWidth',ld)
% grid on
% txt = ['avg%: ',num2str(round(errorilmidflow*100,2))];
% t1 = text(0.7,30,txt);
% t1.FontSize = fontxt;
% txt = ['max%: ',num2str(round(emaxilmidflow*100,2))];
% t1 = text(0.7,23.5,txt);
% t1.FontSize = fontxt;
% txt = ['sys%: ',num2str(round(esysilmidflow*100,2))];
% t1 = text(0.7,17,txt);
% t1.FontSize = fontxt;
% txt = ['dias%: ',num2str(round(ediasilmidflow*100,2))];
% t1 = text(0.7,10.5,txt);
% t1.FontSize = fontxt;
% t2 = legend('in-house','3-D');
% t2.FontSize = fontxt2;
% t2 = title('Volume flow rate - middle of the iliac');
% t2.FontSize = fontxt2;
% t2 = xlabel('time(s)');
% t2.FontSize = fontxt2;
% t2 = ylabel('Q(ml/s)');
% t2.FontSize = fontxt2;
% axis([0 1.2,-20 50]);
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 0.5;
% ax.Layer = 'top';
% ax.XAxis.FontSize = fontxt2;
% ax.YAxis.FontSize = fontxt2;

figure
plot(t,q2ilout*10^6,colour1,'LineWidth',ld)
hold on
plot(outlet_il_flow(1:end,1),outlet_il_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
txt = ['avg%: ',num2str(round(erroriloutflow*100,2))];
t1 = text(0.7,30,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutflow*100,2))];
t1 = text(0.7,23.5,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutflow*100,2))];
t1 = text(0.7,17,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutflow*100,2))];
t1 = text(0.7,10.5,txt);
t1.FontSize = fontxt;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Volume flow rate - outlet of the iliac');
t2.FontSize = fontxt2;
t2 = xlabel('time(s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q(ml/s)');
t2.FontSize = fontxt2;
axis([0 1.2,-20 50]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

%Pressure
figure
plot(t,p2ilout*10^-3,colour1,'LineWidth',ld)
hold on
plot(outlet_il_Pressure(1:end,1),outlet_il_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - outlet of the iliac');
t2.FontSize = fontxt2;
t2 = xlabel('time(s)');
t2.FontSize = fontxt2;
t2 = ylabel('P(kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(erroriloutpressure*100,2))];
t1 = text(0.7,18,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutpressure*100,2))];
t1 = text(0.7,16.8,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutpressure*100,2))];
t1 = text(0.7,15.6,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutpressure*100,2))];
t1 = text(0.7,14.4,txt);
t1.FontSize = fontxt;
grid on
axis([0 1.2,8.5 22]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

% figure
% plot(t,p2mid*10^-3,colour1,'LineWidth',ld)
% hold on
% plot(mid_il_Pressure(1:end,1),mid_il_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
% txt = ['avg%: ',num2str(round(errorilmidpressure*100,2))];
% t1 = text(0.7,18,txt);
% t1.FontSize = fontxt;
% txt = ['max%: ',num2str(round(emaxilmidpressure*100,2))];
% t1 = text(0.7,16.8,txt);
% t1.FontSize = fontxt;
% txt = ['sys%: ',num2str(round(esysilmidpressure*100,2))];
% t1 = text(0.7,15.6,txt);
% t1.FontSize = fontxt;
% txt = ['dias%: ',num2str(round(ediasilmidpressure*100,2))];
% t1 = text(0.7,14.4,txt);
% t1.FontSize = fontxt;
% t2 = legend('in-house','3-D');
% t2.FontSize = fontxt2;
% t2 = title('Pressure - middle of the iliac');
% t2.FontSize = fontxt2;
% t2 = xlabel('time(s)');
% t2.FontSize = fontxt2;
% t2 = ylabel('P(kPa)');
% t2.FontSize = fontxt2;
% grid on
% axis([0 1.2,8.5 22]);
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 0.5;
% ax.Layer = 'top';
% ax.XAxis.FontSize = fontxt2;
% ax.YAxis.FontSize = fontxt2;


figure
plot(t,p1junc*10^-3,colour1,'LineWidth',ld)
hold on
plot(junction_Pressure(1:end,1),junction_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
box on
ax = gca;
ax.LineWidth = 2.5;
txt = ['avg%: ',num2str(round(errorjuncpressure*100,2))];
t1 = text(0.7,18,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxjuncpressure*100,2))];
t1 = text(0.7,16.8,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysjuncpressure*100,2))];
t1 = text(0.7,15.6,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasjuncpressure*100,2))];
t1 = text(0.7,14.4,txt);
t1.FontSize = fontxt;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - junction');
t2.FontSize = fontxt2;
t2 = xlabel('time(s)');
t2.FontSize = fontxt2;
t2 = ylabel('P(kPa)');
t2.FontSize = fontxt2;
grid on
axis([0 1.2,8.5 22]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;

% figure
% plot(t,p1mid*10^-3,colour1,'LineWidth',ld)
% hold on
% plot(mid_ao_Pressure(1:end,1),mid_ao_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
% txt = ['avg%: ',num2str(round(erroraomidpressure*100,2))];
% t1 = text(0.7,18,txt);
% t1.FontSize = fontxt;
% txt = ['max%: ',num2str(round(emaxaomidpressure*100,2))];
% t1 = text(0.7,16.8,txt);
% t1.FontSize = fontxt;
% txt = ['sys%: ',num2str(round(esysaomidpressure*100,2))];
% t1 = text(0.7,15.6,txt);
% t1.FontSize = fontxt;
% txt = ['dias%: ',num2str(round(ediasaomidpressure*100,2))];
% t1 = text(0.7,14.4,txt);
% t1.FontSize = fontxt;
% t2 = legend('in-house','3-D');
% t2.FontSize = fontxt2;
% t2 = title('Pressure - middle of the aorta');
% t2.FontSize = fontxt2;
% t2 = xlabel('time(s)');
% t2.FontSize = fontxt2;
% t2 = ylabel('P(kPa)');
% t2.FontSize = fontxt2;
% grid on
% axis([0 1.2,8.5 22]);
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 0.5;
% ax.Layer = 'top';
% ax.XAxis.FontSize = fontxt2;
% ax.YAxis.FontSize = fontxt2;

figure
plot(t,p1in*10^-3,colour1,'LineWidth',ld)
hold on
plot(Inlet_ao_Pressure(1:end,1),Inlet_ao_Pressure(1:end,2)*10^-3,colour2,'LineWidth',ld);
box on
ax = gca;
ax.LineWidth = 2.5;
txt = ['avg%: ',num2str(round(errorinletpressure*100,2))];
t1 = text(0.7,18,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxinletpressure*100,2))];
t1 = text(0.7,16.8,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysinletpressure*100,2))];
t1 = text(0.7,15.6,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasinletpressure*100,2))];
t1 = text(0.7,14.4,txt);
t1.FontSize = fontxt;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt2;
t2 = title('Pressure - inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time(s)');
t2.FontSize = fontxt2;
t2 = ylabel('P(kPa)');
t2.FontSize = fontxt2;
grid on
axis([0 1.2,8.5 22]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;