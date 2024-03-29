%% This is programme to calculate for the common carrotid artery
clc;
clear;

%% Defining the constant parameters of the problem
a = 0.001;                                          %Very small angle

%Properties
L = 126*10^-3;                                      %Table 1 (Flores 2016)
R = 3*10^-3;                                        %Table 1 (Flores 2016)
h = 0.3*10^-3;                                      %Table 1 (Flores 2016)
rho = 1060;                                         %Table 1 (Flores 2016)
E = 700*10^3;                                       %Table 1 (Flores 2016)
v = 0.5;        
RW1 = 2.4875*10^8;                                  %Table 1 (Flores 2016)
RW2 = 1.8697*10^9;                                  %Table 1 (Flores 2016)
Cwk = 1.7529*10^-10;                                %Table 1 (Flores 2016)

% L = 0.139467775850962;
% R = 0.00298937029503352;
% bet = 0.00341884603442744;
%                                      %Table 1 (Flores 2016)
% rho = 1060;                                         %Table 1 (Flores 2016)
% E = 700*10^3;                                       %Table 1 (Flores 2016)
% v = 0.5;
% h = (1-v^2)/(E*bet); 
% RW1 = 273076051.267355;                                  %Table 1 (Flores 2016)
% RW2 = 1852305574.71689;                                  %Table 1 (Flores 2016)
% Cwk = 1.74814291325816e-10;                                %Table 1 (Flores 2016)

%% Importing the data from the Flores plots
load('automeris/common_carotid_artery/CC_inlet_BC.csv')
load('automeris/common_carotid_artery/inlet_pressure.csv')
load('automeris/common_carotid_artery/mid_pressure.csv')
load('automeris/common_carotid_artery/outlet_pressure.csv')
load('automeris/common_carotid_artery/pin_pout.csv')
load('grabit/CC_mid.mat')
load('grabit/CC_outlet.mat')

x = CC_inlet_BC(1:end,1);                            %Extracting the time data to a vector
y = CC_inlet_BC(1:end,2)*10^-6;                      %Extracting the Q data to a vector
N = length(y);
T = CC_inlet_BC(end,1);                              %Extracting the period of th signal
xi = (x(1):T/((N)):T)';                              %Creating an equispaced time vector to be used in the fft
yi = interp1q(x,y,xi);                               %Interpolating the values of the Q to match the new time vector
Qin = yi(1:end);                                     %Vector with the values of the BC
t = xi(1:end);
N = length(Qin);

%% Fourier Transform
F = fft(Qin(1:N))/N;

% Define the fundamental harmonic omega and the number of harmonics
om = 2*pi/T;
nh = N/2;

%% Frequency domain calculations
for ih = 0:nh  
     
  if ih == 0 
      
      omega(ih+1) = 10^-10;                               %Assining a very small value to ω to avoid singularity
      
  else
      omega(ih+1) = -ih*om;                               
      
  end
  ih=ih+1;
%% Backward Propagation 
    
    %Impedance Z 
    Z = (RW1+RW2-1i*omega(ih)*RW1*RW2*Cwk)/(1-1i*omega(ih)*RW2*Cwk);

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
  
%% Forward Propagation     
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
     
        
      q1 = q1 + real(Q1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      p1 = p1 + real(P1(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
        
      q1mid = q1mid + real(Q1mid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      p1mid = p1mid + real(P1mid(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      
      q1outlet = q1outlet + real(Q1outlet(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));
      p1outlet = p1outlet + real(P1outlet(ih+1)*2*F(ih+1)*exp(-1i*omega(ih+1)*t));

 end
 

%% Error calculations
%Average

yiinletpressure = interp1q(inlet_pressure(1:end,1),inlet_pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errinlet(i) = abs((yiinletpressure(i) - p1(i))/(yiinletpressure(i)))^2;

end
errorinletpressure = sqrt(mean(errinlet));

yimidpressure = interp1q(mid_pressure(1:end,1),mid_pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errmid(i) = abs((yimidpressure(i) - p1mid(i))/(yimidpressure(i)))^2;

end
errormidpressure = sqrt(mean(errmid));

yimidflow = interp1q(CC_mid(1:end,1),CC_mid(1:end,2),xi);
for i = 2:length(t)-2
    
     errmid(i) = abs((yimidflow(i)*10^-6 - q1mid(i))/(max(yimidflow)*10^-6))^2;

end
errormidflow = sqrt(mean(errmid));

yioutletpressure = interp1q(outlet_pressure(1:end,1),outlet_pressure(1:end,2),xi);
for i = 2:length(t)-2
    
    errmid(i) = abs((yioutletpressure(i) - p1outlet(i))/(yioutletpressure(i)))^2;

end
erroroutletpressure = sqrt(mean(errmid));

yioutletflow = interp1q(CC_outlet(1:end,1),CC_outlet(1:end,2),xi);
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

% Ploting the results
ld = 2.5;
colour1 = 'b';
colour2 = 'r';
fontxt = 21;
fontxt2 = 21;
fontxt3 = 21;

figure(1)

plot(t,Qin*10^6,colour1,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = title('Volume flow rate - Inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
axis([0 1.2,2 14]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;

figure
plot(t,q1mid*10^6,colour1,'LineWidth',ld)
hold on
plot(CC_mid(1:end,1),CC_mid(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt3;
t2 = title('Volume flow rate - Midpoint');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(errormidflow*100,2))];
t1 = text(0.8,10.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxmidflow*100,2))];
t1 = text(0.8,9.5,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysmidflow*100,2))];
t1 = text(0.8,8.5,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasmidflow*100,2))];
t1 = text(0.8,7.5,txt);
t1.FontSize = fontxt;
axis([0 1.2,2 14]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;

figure
plot(t,q1outlet*10^6,colour1,'LineWidth',ld)
hold on
plot(CC_outlet(1:end,1),CC_outlet(1:end,2),colour2,'LineWidth',ld)
grid on
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt3;
t2 = title('Volume flow rate - Outlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q (ml/s)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(erroroutletflow*100,2))];
t1 = text(0.8,10.5,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutflow*100,2))];
t1 = text(0.8,9.5,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutflow*100,2))];
t1 = text(0.8,8.5,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutflow*100,2))];
t1 = text(0.8,7.5,txt);
t1.FontSize = fontxt;
axis([0 1.2,2 14]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;

figure
plot(t,p1*10^-3,colour1,'LineWidth',ld)
hold on
box on
ax = gca;
ax.LineWidth = 2.5;
plot(inlet_pressure(1:end,1),inlet_pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt3;
t2 = title('Pressure - Inlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(errorinletpressure*100,2))];
t1 = text(0.8,16.4,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxinletpressure*100,2))];
t1 = text(0.8,15.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysinletpressure*100,2))];
t1 = text(0.8,15,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasinletpressure*100,2))];
t1 = text(0.8,14.3,txt);
t1.FontSize = fontxt;
axis([0 1.2,10 19]);
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
box on
ax = gca;
ax.LineWidth = 2.5;
plot(mid_pressure(1:end,1),mid_pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt3;
t2 = title('Pressure - Midpoint');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(errormidpressure*100,2))];
t1 = text(0.8,16.4,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxmidpressure*100,2))];
t1 = text(0.8,15.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysmidpressure*100,2))];
t1 = text(0.8,15,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasmidpressure*100,2))];
t1 = text(0.8,14.3,txt);
t1.FontSize = fontxt;
axis([0 1.2,10 19]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;
grid on

figure
plot(t,p1outlet*10^-3,colour1,'LineWidth',ld)
hold on
plot(outlet_pressure(1:end,1),outlet_pressure(1:end,2)*10^-3,colour2,'LineWidth',ld)
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt3;
t2 = title('Pressure - Outlet');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(erroroutletpressure*100,2))];
t1 = text(0.8,16.4,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutpressure*100,2))];
t1 = text(0.8,15.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutpressure*100,2))];
t1 = text(0.8,15,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutpressure*100,2))];
t1 = text(0.8,14.3,txt);
t1.FontSize = fontxt;
axis([0 1.2,10 19]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;
grid on


figure
plot(t,(p1-p1outlet)*10^-3,colour1,'LineWidth',ld)
hold on
plot(pin_pout(1:end,1),pin_pout(1:end,2)*10^-3,colour2,'LineWidth',ld)
box on
ax = gca;
ax.LineWidth = 2.5;
t2 = legend('1-D (Present)','3-D');
t2.FontSize = fontxt3;
t2 = title('Pin - Pout');
t2.FontSize = fontxt2;
t2 = xlabel('time (s)');
t2.FontSize = fontxt2;
t2 = ylabel('P (kPa)');
t2.FontSize = fontxt2;
txt = ['avg%: ',num2str(round(erroroutletpressure*100,2))];
t1 = text(0.8,16.4,txt);
t1.FontSize = fontxt;
txt = ['max%: ',num2str(round(emaxoutpressure*100,2))];
t1 = text(0.8,15.7,txt);
t1.FontSize = fontxt;
txt = ['sys%: ',num2str(round(esysoutpressure*100,2))];
t1 = text(0.8,15,txt);
t1.FontSize = fontxt;
txt = ['dias%: ',num2str(round(ediasoutpressure*100,2))];
t1 = text(0.8,14.3,txt);
t1.FontSize = fontxt;
axis([0 1.2,-0.4 1]);
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
ax.XAxis.FontSize = fontxt3;
ax.YAxis.FontSize = fontxt3;
grid on