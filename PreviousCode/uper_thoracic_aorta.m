%% This is programme to calculate for the uper thoracic aorta
clc;
clear;

%% Defining the constant parameters of the problem
a = 0.001;  

%Properties
L = 24.137*10^-2;                                   %Table 2 (Flores 2016)
R = 1.27*10^-2;                                     %Table 2 (Flores 2016)
h = 1.2*10^-3;                                      %Table 2 (Flores 2016)
rho = 1060;                                         %Table 2 (Flores 2016)
E = 400*10^3;                                       %Table 2 (Flores 2016)
v = 0.5;        
RW1 = 1.1752*10^7;                                  %Table 2 (Flores 2016)
RW2 = 1.1167*10^8;                                  %Table 2 (Flores 2016)
Cwk = 1.0163*10^-8;                                 %Table 2 (Flores 2016)

% L = 0.231745693419459;
% R = 0.0117542467502781;
% RW1 = 11672188.9876804;
% RW2 = 110825572.563458;
% Cwk = 1.01406378445226e-08;
% bet = 0.00174707638017433;
% E = 400*10^3;
% v = 0.5;        
% h = (1-v^2)/(E*bet);

%% Importing the data from the Flores plots
load('automeris/uper_thoracic_aorta/UTA_BC.csv')
load('automeris/uper_thoracic_aorta/Inlet_pressure.csv')
load('automeris/uper_thoracic_aorta/midpoint_pressure.csv')
load('automeris/uper_thoracic_aorta/outlet_pressure.csv')
load('automeris/uper_thoracic_aorta/midpoint_flow.csv')
load('automeris/uper_thoracic_aorta/outlet_flow.csv')

x = UTA_BC(1:end,1);
y = UTA_BC(1:end,2)*10^-6;
N = length(y);
T = UTA_BC(end,1);
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
      omega(ih+1) = ih*om;
      
  end
  ih=ih+1;
    %% Step 1
    %Firstly the term B1/A1 
    
    %Impedance Z 
    Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);
    %Z = -real(Z1)+imag(Z1)*1i;
    %Bessel functions at s1out
    s1out = (R-tan(a)*L)/sin(a);
    [J_13_s1out,Y_13s1out,J_43s1out,Y_43s1out,fs1out] = besselfunctions(a,s1out,omega(ih),E,rho,v,h); 
    Y_s1out = (2*pi*(1-cos(a)))*(fs1out/rho)^0.5*s1out^2.5;
      
    
    B1_A1 = -(J_13_s1out+1i*Y_s1out*J_43s1out*Z)/(Y_13s1out+1i*Y_s1out*Y_43s1out*Z);
    
    s1in = (R-tan(a)*0)/sin(a);
    [J_13_s1in,Y_13s1in,J_43s1in,Y_43s1in,fs1in] = besselfunctions(a,s1in,omega(ih),E,rho,v,h); 
    Y_s1in = (2*pi*(1-cos(a)))*(fs1in/rho)^0.5*s1in^2.5;
    A1tilda = -1/(1i*Y_s1in*(s1in^-0.5)*(J_43s1in+B1_A1*Y_43s1in));
    B1tilda = B1_A1*A1tilda;
    
    x1 = 0;
    s1 = (R-tan(a)*x1)/sin(a);
    [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E,rho,v,h); 
    Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
    Q1(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
    P1(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
    
    x1 = L/2;
    s1 = (R-tan(a)*x1)/sin(a);
    [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E,rho,v,h); 
    Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
    Q1mid(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
    P1mid(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
    
    x1 = L;
    s1 = (R-tan(a)*x1)/sin(a);
    [J_13s1,Y_13s1,J_43s1,Y_43s1,fs1] = besselfunctions(a,s1,omega(ih),E,rho,v,h); 
    Y_s1 = (2*pi*(1-cos(a)))*(fs1/rho)^0.5*s1^2.5;
    Q1outlet(ih) = -(1i*Y_s1*(s1^-0.5)*(A1tilda*J_43s1+B1tilda*Y_43s1));
    P1outlet(ih) = ((s1^-0.5)*(A1tilda*J_13s1+B1tilda*Y_13s1));
    
end

%% Inverse Fourier
q1 = zeros(size(t)) + real(Q1(1)*F(1));
p1 = zeros(size(t)) + real(P1(1)*F(1));
q1mid = zeros(size(t)) + real(Q1mid(1)*F(1));
p1mid = zeros(size(t)) + real(P1mid(1)*F(1));
q1outlet = zeros(size(t)) + real(Q1outlet(1)*F(1));
p1outlet = zeros(size(t)) + real(P1outlet(1)*F(1));

 for ih = 1:nh
     
        
      q1 = q1 + real(conj(Q1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p1 = p1 + real(conj(P1(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
        
      q1mid = q1mid + real(conj(Q1mid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p1mid = p1mid + real(conj(P1mid(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      
      q1outlet = q1outlet + real(conj(Q1outlet(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));
      p1outlet = p1outlet + real(conj(P1outlet(ih+1))*2*F(ih+1)*exp(1i*omega(ih+1)*t));

 end
 

 
%% Error calculations
%Average
yiinletpressure = interp1q(Inlet_pressure(1:end,1),Inlet_pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errinlet(i) = abs((yiinletpressure(i) - p1(i))/(yiinletpressure(i)))^2;

end
errorinletpressure = sqrt(mean(errinlet));

yimidpressure = interp1q(midpoint_pressure(1:end,1),midpoint_pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errmid(i) = abs((yimidpressure(i) - p1mid(i))/(yimidpressure(i)))^2;

end
errormidpressure = sqrt(mean(errmid));

yimidflow = interp1q(midpoint_flow(1:end,1),midpoint_flow(1:end,2),xi);
for i = 2:length(t)-2
    
     errmid(i) = abs((yimidflow(i)*10^-6 - q1mid(i))/(max(yimidflow)*10^-6))^2;

end
errormidflow = sqrt(mean(errmid));

yioutletpressure = interp1q(outlet_pressure(1:end,1),outlet_pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errmid(i) = abs((yioutletpressure(i) - p1outlet(i))/(yioutletpressure(i)))^2;

end
erroroutletpressure = sqrt(mean(errmid));

yioutletflow = interp1q(outlet_flow(1:end,1),outlet_flow(1:end,2),xi);
for i = 2:length(t)-2
    
     errmid(i) = abs((yioutletflow(i)*10^-6 - q1outlet(i))/(max(yioutletflow)*10^-6))^2;

end
erroroutletflow = sqrt(mean(errmid));

%MAX
emaxmidflow = max(abs((q1mid(2:end-1)-yimidflow(2:end-1)*10^-6)/max(yimidflow(2:end-1)*10^-6)));
emaxoutflow = max(abs((q1outlet(2:end-1)-yioutletflow(2:end-1)*10^-6)/max(yioutletflow(2:end-1)*10^-6)));
emaxinletpressure = max(abs((p1(2:end-1)-yiinletpressure(2:end-1))/max(yiinletpressure(2:end-1))));
emaxmidpressure = max(abs((p1mid(2:end-1)-yimidpressure(2:end-1))/max(yimidpressure(2:end-1))));
emaxoutpressure = max(abs((p1outlet(2:end-1)-yioutletpressure(2:end-1))/max(yioutletpressure(2:end-1))));


%SYS
esysmidflow = (max(q1mid(2:end-1)-max(yimidflow(2:end-1))*10^-6))/max(yimidflow(2:end-1)*10^-6);
esysoutflow = (max(q1outlet(2:end-1)-max(yioutletflow(2:end-1))*10^-6))/max(yioutletflow(2:end-1)*10^-6);
esysinletpressure = (max(p1(2:end-1)-max(yiinletpressure(2:end-1))))/max(yiinletpressure(2:end-1));
esysmidpressure = (max(p1mid(2:end-1)-max(yimidpressure(2:end-1))))/max(yimidpressure(2:end-1));
esysoutpressure = (max(p1outlet(2:end-1)-max(yioutletpressure(2:end-1))))/max(yioutletpressure(2:end-1));

%DIAS
ediasmidflow = (min(q1mid(2:end-1)-min(yimidflow(2:end-1))*10^-6))/max(yimidflow(2:end-1)*10^-6);
ediasoutflow = (min(q1outlet(2:end-1)-min(yioutletflow(2:end-1))*10^-6))/max(yioutletflow(2:end-1)*10^-6);
ediasinletpressure = (min(p1(2:end-1)-min(yiinletpressure(2:end-1))))/max(yiinletpressure(2:end-1));
ediasmidpressure = (min(p1mid(2:end-1)-min(yimidpressure(2:end-1))))/max(yimidpressure(2:end-1));
ediasoutpressure = (min(p1outlet(2:end-1)-min(yioutletpressure(2:end-1))))/max(yioutletpressure(2:end-1));



 
%% Ploting the results

currmodel = '1-D (Known)';

ld = 2.5;
colour1 = 'b';
colour2 = 'r';
fontxt = 21;
fontxt2 = 21;
fontxt3 = 21;
figure
plot(t,Qin*10^6,colour1,'LineWidth',ld)
box on
ax = gca;
ax.LineWidth = 2.5;
grid on
t2 = title('Volume flow rate - Inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
axis([0 1,-100 550]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;

figure
plot(t,q1mid*10^6,colour1,'LineWidth',ld)
hold on
plot(midpoint_flow(1:end,1),midpoint_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend(currmodel,'3-D');
t2.FontSize = fontxt3;
t2 = title('Volume flow rate - Midpoint');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(errormidflow*100,2))];
t1 = text(0.67,350,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxmidflow*100,2))];
t1 = text(0.67,300,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysmidflow*100,2))];
t1 = text(0.67,250,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasmidflow*100,2))];
t1 = text(0.67,200,txt);
t1.FontSize = fontxt;
axis([0 1,-100 550]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;

figure
plot(t,q1outlet*10^6,colour1,'LineWidth',ld)
hold on
plot(outlet_flow(1:end,1),outlet_flow(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend(currmodel,'3-D');
t2.FontSize = fontxt3;
t2 = title('Volume flow rate - Outlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(erroroutletflow*100,2))];
t1 = text(0.67,350,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutflow*100,2))];
t1 = text(0.67,300,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutflow*100,2))];
t1 = text(0.67,250,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutflow*100,2))];
t1 = text(0.67,200,txt);
t1.FontSize = fontxt;
axis([0 1,-100 550]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;
% 
figure
plot(t,p1*10^-3,colour1,'LineWidth',ld)
hold on
plot(Inlet_pressure(1:end,1),Inlet_pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend(currmodel,'3-D');
t2.FontSize = fontxt3;
t2 = title('Pressure - Inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(errorinletpressure*100,2))];
t1 = text(0.67,15.7,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxinletpressure*100,2))];
t1 = text(0.67,14.8,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysinletpressure*100,2))];
t1 = text(0.67,13.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasinletpressure*100,2))];
t1 = text(0.67,13,txt);
t1.FontSize = fontxt;
axis([0 1,8 19]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;
grid on

figure
plot(t,p1mid*10^-3,colour1,'LineWidth',ld)
hold on
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
plot(midpoint_pressure(1:end,1),midpoint_pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
t2 = legend(currmodel,'3-D');
t2.FontSize = fontxt3;
t2 = title('Pressure - Midpoint');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(errormidpressure*100,2))];
t1 = text(0.67,15.7,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxmidpressure*100,2))];
t1 = text(0.67,14.8,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysmidpressure*100,2))];
t1 = text(0.67,13.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasmidpressure*100,2))];
t1 = text(0.67,13,txt);
t1.FontSize = fontxt;
axis([0 1,8 19]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;


figure
plot(t,p1outlet*10^-3,colour1,'LineWidth',ld)
hold on
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
plot(outlet_pressure(1:end,1),outlet_pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
t2 = legend(currmodel,'3-D');
t2.FontSize = fontxt3;
t2 = title('Pressure - Outlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(erroroutletpressure*100,2))];
t1 = text(0.67,15.7,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutpressure*100,2))];
t1 = text(0.67,14.8,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutpressure*100,2))];
t1 = text(0.67,13.9,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutpressure*100,2))];
t1 = text(0.67,13,txt);
t1.FontSize = fontxt;
axis([0 1,8 19]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;

