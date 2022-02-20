%% Programme for the Flores implementation of the model bifurcation case
clear;
clc;
%% Defining the constant parameters of the problem

%Parent vessel (1)
E1 = 500*10^3;                                       %Table 3 (Flores 2016)
rho = 1060;                                         %Table 3 (Flores 2016)
v1 = 0.5;                                            %Table 3 (Flores 2016)
h1 = 1.032*10^-3;                                    %Table 3 (Flores 2016)
R1 = 0.89*10^-2;
L1 = 8.6*10^-2;
%Daughter vessels (2) and (3)
E2 = 700*10^3;                                       %Table 3 (Flores 2016)
rho = 1060;                                         %Table 3 (Flores 2016)
v2 = 0.5;                                            %Table 3 (Flores 2016)
h2 = 0.72*10^-3;                                     %Table 3 (Flores 2016)
R2 = 0.6125*10^-2;
L2 = 8.5*10^-2;


%Impedence calculation data
RW1 = 6.8123*10^7;
RW2 = 3.1013*10^9;
Cwk = 3.6664*10^-10;

%Importing the data from the Flores plots
load('automeris/bifurcation/BC.csv')
load('automeris/bifurcation/Inlet_ao_Pressure.csv')
load('automeris/bifurcation/mid_ao_Pressure.csv')
load('grabit/ao_mid_flow.mat')
load('automeris/bifurcation/flores/junc_flow.csv')
load('grabit/il_mid_flow.mat')
load('grabit/il_outlet_flow.mat')


x = BC(1:end,1);
y = BC(1:end,2)*10^-6;
xmid = ao_mid_flow(1:end,1);
ymid = ao_mid_flow(1:end,2)*10^-6;
% xjunc = junction_flow(1:end,1);
% yjunc = junction_flow(1:end,2);
xilmid = il_mid_flow(1:end,1);
yilmid = il_mid_flow(1:end,2)*10^-6;
xiout = il_outlet_flow(1:end,1);
yiout = il_outlet_flow(1:end,2)*10^-6;
N = length(y);
T = BC(end,1);
xi = (x(1):T/((2*N)):T)';
yi = interp1q(x,y,xi);
Qin = yi(1:end);
t = xi(1:end);
N = length(Qin);
F = fft(Qin(1:N))/N;


% Define the fundamental harmonic omega and the number of
% harmonics
om = 2*pi/T;
nh = 20;

%Initialization of vectors
for ih = 0:nh   %Work with the first 10 harmonics as in Papadakis 2019
%% Backward Propagation       
     
    if ih == 0
        omega(ih+1) = 10^-10;
    
    else
         omega(ih+1) = ih*om;
    end
    ih=ih+1;
   
    Z = (RW1+RW2-1i*(omega(ih))*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk); 
    R01 = R1;
    Rd1 = 0.86*10^-2;
    A01 = pi*R01^2;
    hta1 = 4*10^-3;
    C1 = (3*pi*R01*Rd1^2)/(2*E1*h1);
    k1 = sqrt(rho*1i*omega(ih)/hta1);
    K_omega1 = -(hta1/(1i*omega(ih)*rho))*(1-(2*besselj(1,(k1*R01)))/(k1*R01*besselj(0,(k1*R01))));
    M1 = sqrt(1i*omega(ih)*C1*(A01)*K_omega1/hta1);
    Kc1 = sqrt(1i*omega(ih)*C1*hta1/(A01*K_omega1));
    
    R02 = R2;
    Rd2 = 0.6*10^-2;
    A02 = pi*R02^2;
    hta2 = 4*10^-3;
    C2 = (3*pi*R02*Rd2^2)/(2*E2*h2);
    k2 = sqrt(rho*1i*omega(ih)/hta2);
    K_omega2 = -(hta2/(1i*omega(ih)*rho))*(1-(2*besselj(1,(k2*R02)))/(k2*R02*besselj(0,(k2*R02))));
    M2 = sqrt(1i*omega(ih)*C2*(A02)*K_omega2/hta2);
    Kc2 = sqrt(1i*omega(ih)*C2*hta2/(A02*K_omega2));
    
    x = 0;
    k12 = M2*cos(Kc2*L2)/sin(Kc2*L2);
    k31 = M1*sin(Kc1*L1)/cos(Kc1*L1);
    k22 = M2*(1/sin(Kc2*L2));
    P1fl(ih) = -((1+Z*k12)*cos(Kc1*x))/(cos(Kc1*L1)^2*(k31-2*k12+Z*(k12*k31-2*k12^2+2*k22^2)))+(sin(Kc1*L1)/(M1*cos(Kc1*L1)))*cos(Kc1*x)-sin(Kc1*x)/M1;   
    Q1fl(ih) = -((1+Z*k12)*sin(Kc1*x)*M1)/(cos(Kc1*L1)^2*(k31-2*k12+Z*(k12*k31-2*k12^2+2*k22^2)))+(sin(Kc1*L1)/(cos(Kc1*L1)))*sin(Kc1*x)+cos(Kc1*x);
    
    x2 = L2;
    P2fl(ih) = -((1+Z*k12)*(sin(Kc2*L2)*cos(Kc2*x2)-cos(Kc2*L2)*sin(Kc2*x2))+Z*k22*sin(Kc2*x2))/(sin(Kc2*L2)*cos(Kc2*L2)*(k31-2*k12+Z*(k12*k31-2*k12^2+2*k22^2)));
    Q2fl(ih) = -M2*((1+Z*k12)*(sin(Kc2*L2)*sin(Kc2*x2)+cos(Kc2*L2)*cos(Kc2*x2))-Z*k22*cos(Kc2*x2))/(sin(Kc2*L2)*cos(Kc2*L2)*(k31-2*k12+Z*(k12*k31-2*k12^2+2*k22^2)));

end

p1fl = zeros(size(t)) + real(conj(P1fl(1))*F(1));
q1fl = zeros(size(t)) + real(conj(Q1fl(1))*F(1));
p2fl = zeros(size(t)) + real(conj(P2fl(1))*F(1));
q2fl = zeros(size(t)) + real(conj(Q2fl(1))*F(1));



 for ih = 1:nh
     
     p1fl = p1fl + real(conj(P1fl(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
     q1fl = q1fl+ real(conj(Q1fl(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
     p2fl = p2fl + real(conj(P2fl(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
     q2fl = q2fl+ real(conj(Q2fl(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));

 end

p1fl = q2fl;
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

%Impedence calculation data
RW1 = 6.8123*10^7;
RW2 = 3.1013*10^9;
Cwk = 3.6664*10^-10;


%% Frequency domain calculations
for ih = 0:nh   
%% Backward Propagation       
     
    if ih == 0
        omega(ih+1) = 10^-10;
    
    else
         omega(ih+1) = ih*om;
    end
    
    ih=ih+1;
    
    %% Step 1
    %Firstly the terms B2/A2 and B3/A3 need to be calculated
    
    %Impedance Z (same for both daughter vessels)
    Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);
    %Bessel functions at s2out and s3out
    s2out = (R2-tan(a2)*L2)/sin(a2);
    s3out = s2out;
    [J_13_s2out,Y_13s2out,J_43s2out,Y_43s2out,fs2out] = besselfunctions(a2,s2out,omega(ih),E2,rho2,v2,h2); 
    J_13_s3out = J_13_s2out;
    Y_13s3out = Y_13s2out;
    J_43s3out = J_43s2out;
    Y_43s3out = Y_43s2out;
    fs3out = fs2out;
    
    Y_s2out = (2*pi*(1-cos(a2)))*(fs2out/rho2)^0.5*s2out^2.5;
    
    B2_A2 = -(J_13_s2out+1i*Y_s2out*J_43s2out*Z)/(Y_13s2out+1i*Y_s2out*Y_43s2out*Z);
    
    %% Step 2
    %Calculation of Yeff(s2in) and Yeff(s3in)
    s2in = R2/sin(a2);
    [J_13_s2in,Y_13s2in,J_43s2in,Y_43s2in,fs2in] = besselfunctions(a2,s2in,omega(ih),E2,rho2,v2,h2); 
    J_13_s3in = J_13_s2in;
    Y_13s3in = Y_13s2in;
    J_43s3in = J_43s2in;
    Y_43s3in = Y_43s2in;
    Y_s2in = (2*pi*(1-cos(a2)))*(fs2in/rho2)^0.5*s2in^2.5;
    
    Yeff_s2in = -1i*Y_s2in*(J_43s2in+B2_A2*Y_43s2in)/(J_13_s2in+B2_A2*Y_13s2in);
    Yeff_s3in = Yeff_s2in;
    
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
    %Q1junc(ih) = (Y_s1out*(s1out^-0.5)*(J_43s1out+B1_A1*Y_43s1out))/(Y_s0*(s0^(-0.5))*(J_43s0+B1_A1*Y_43s0));
    Q1junc(ih) = -1i*(Y_s1out*(s1out^-0.5))*(A1tilda*J_43s1out+B1tilda*Y_43s1out);
%     P1junc(ih) = -(s1out^-0.5)*((J_13_s1out)+B1_A1*Y_13s1out)/(1i*Y_s0*(s0^(-0.5)*(J_43s0+B1_A1*Y_43s0)));

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
%    Q2tot(ih) = 2*P2in(ih)*Yeff_s2in;
%    P2in2(ih) = (Q1junc(ih)/(Yeff_s3in+Yeff_s2in)); 
   
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
     
       q1 = q1 + real(Q1(ih+1)*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       q1mid = q1mid + real(conj(Q1mid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       q1junc = q1junc + real(conj(Q1junc(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       q2ilmid = q2ilmid + real(conj(Q2ilmid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       q2ilout = q2ilout + real(conj(Q2ilout(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       p1in = p1in + real(conj(P1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       p1mid = p1mid + real(conj(P1mid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       p1junc = p1junc + real(conj(P1junc(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       p2in = p2in + real(conj(P2in(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       p2mid = p2mid + real(conj(P2ilmid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
       p2ilout = p2ilout + real(conj(P2ilout(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
   
 end

 
% p1in = p1in;
figure

plot(t,q2ilout*10^6,'b','LineWidth',2.5)
% plot(t,Qin*10^6,'b','LineWidth',2.5)
hold on
plot(t,q2fl*10^6,'--k','LineWidth',2.5)



for i = 1:length(t)
    
     erriloutf(i) = abs((q2fl(i) - q2ilout(i))/(max(q2fl)))^2;

end
erroriloutflow = sqrt(mean(erriloutf));

emaxoutflow = max(abs((q2fl-q2ilout)/max(q2fl)));

esysoutflow = (max(q2fl-max(q2ilout)))/max(q2fl);

ediasoutflow = (min(q2fl-min(q2ilout)))/max(q2fl);

% %AVG
% for i = 1:length(t)
%     
%     errinletp(i) = abs((p1fl(i) - p1in(i))/(p1fl(i)))^2;
% 
% end
% 
% errorinletpressure = sqrt(mean(errinletp));
% 
% 
% %MAX
% emaxinletpressure = max(abs((p1fl-p1in)/max(p1fl)));
% 
% %SYS
% esysinletpressure = (max(p1fl)-max(p1in))/max(p1fl);
% 
% %DIAS
% ediasinletpressure = (min(p1fl-min(p1in)))/max(p1fl);


fontxt = 21;
fontxt2 = 21;

grid on
box on
ax = gca;
ax.LineWidth = 2.5;
% t2 = legend('in-house','GDEM');
% t2.FontSize = fontxt2;
t2 = title('Volume flow rate - inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
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
axis([0 1.2,-40 100]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;
