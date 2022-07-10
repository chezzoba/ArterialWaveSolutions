%% This is programme to calculate for the common carrotid artery
clc;
clear;
addpath('../VesselModels/');

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
be = (1 - v^2)/(E*h);

% L = 0.139467775850962;
% R = 0.00298937029503352;
% be = 0.00341884603442744;
%                                      %Table 1 (Flores 2016)
% rho = 1060;                                         %Table 1 (Flores 2016)
% RW1 = 273076051.267355;                                  %Table 1 (Flores 2016)
% RW2 = 1852305574.71689;                                  %Table 1 (Flores 2016)
% Cwk = 1.74814291325816e-10;                                %Table 1 (Flores 2016)

%% Importing the data from the Flores plots
load('../PreviousCode/automeris/common_carotid_artery/CC_inlet_BC.csv')
load('../PreviousCode/automeris/common_carotid_artery/inlet_pressure.csv')
load('../PreviousCode/automeris/common_carotid_artery/mid_pressure.csv')
load('../PreviousCode/automeris/common_carotid_artery/outlet_pressure.csv')
load('../PreviousCode/automeris/common_carotid_artery/pin_pout.csv')
load('../PreviousCode/grabit/CC_mid.mat')
load('../PreviousCode/grabit/CC_outlet.mat')

[omegas, F, t, Qin] = Vessel.ProcessBC(CC_inlet_BC);

artery = Vessel(R, L, a, be, rho, [RW1, RW2, Cwk]);
artery.type = 1;

artery = artery.backpropagate(omegas);
xin = 0;
[Q1, P1] = artery.forwardpropagate(omegas, artery.s(xin));

xmid = artery.L/2;
[Q1mid, P1mid] = artery.forwardpropagate(omegas, artery.s(xmid));

xout = artery.L;
[Q1outlet, P1outlet] = artery.forwardpropagate(omegas, artery.s(xout));

sols = [Q1; P1; Q1mid; P1mid; Q1outlet; P1outlet];
solt = Vessel.InverseFourierTransform(t, omegas, sols, F);

%% Inverse Fourier
q1 = solt(1, :);
p1 = solt(2, :);
q1mid = solt(3, :);
p1mid = solt(4, :);
q1outlet = solt(5, :);
p1outlet = solt(6, :);

xi = t;
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
