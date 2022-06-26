%% Programme for implementation of the aortic arch model
clc;
clear;
addpath('../VesselModels/');

plotting = true;

%% Importing the data from the Flores plots
load('../PreviousCode/automeris/aorta/BC.csv')
load('../PreviousCode/automeris/aorta/Inlet_Pressure.csv')
load('../PreviousCode/automeris/aorta/bc_outlet_11_flow.csv')
load('../PreviousCode/automeris/aorta/bc_outlet_11_Pressure.csv')
load('../PreviousCode/automeris/aorta/lcca_outlet_12_flow.csv')
load('../PreviousCode/automeris/aorta/lcca_outlet_12_Pressure.csv')
load('../PreviousCode/automeris/aorta/lsub_outlet_13_flow.csv')
load('../PreviousCode/automeris/aorta/lsub_outlet_13_Pressure.csv')

rho = 1060; % kg/m^3
seg1 = 0.5; % L1
WKP7 = [59149000, 1.0174e+09, 9.0686e-10];

ves1 = Vessel(0.0152, 0.0704*seg1, 0.0184750925619521, 0.001325687943664,...
    1060, [0, 0, 0]);
ves1.type = 5;

R2 = ves1.R - tan(ves1.a)*seg1*ves1.L;

bi2_8_3 = Bifurcation([R2, 0.00635, 0.0139], [0.0704*(1-seg1), 0.0340, 0.008],...
    [0.00132568794366357, 0.00192990582059595, 0.00140439444384107],...
    rho, [0.0184750925619521, 0.001, 0.02499479361892],...
    51918000, 1060800000, 8.6974e-10);

bi3_9_4 = Bifurcation([0.0139 0.0036 0.0137], [0.0080 0.0340 0.0090],...
    [0.001404394443841,0.002421354408802,0.001412397459944], rho,...
    [0.024994793618920,1e-03,0.022218565326719],...
    191515000, 5.22129e+09, 1.767e-10);

bi4_10_5 = Bifurcation([0.0137,0.0048,0.0135], [0.009,0.034,0.064737],...
    [0.001412397459944,0.002158149171271,0.001388888888889], rho,...
    [0.022218565326719,1e-03,0.018534417520118],...
    98820000, 1.30183e+09, 7.0871e-10);

seg6 = 0.5;
ves6 = Vessel(0.0123, seg6*0.152, 0.015788161735825, 0.001392773178531,...
    rho, [0, 0, 0]);

R7 = ves6.R - tan(ves6.a)*seg6*ves6.L;

ves7 = Vessel(R7, 0.152*(1-seg6), 0.015788161735825, 0.001392773178531,...
    rho, WKP7);
ves7.type = 2;

[omegas, F, t] = Vessel.ProcessBC(BC);

%% Backward Propagation

ves7 = ves7.backpropagate(omegas);

ves6 = ves6.backpropagate(omegas, ves7);

ves5 = bi4_10_5.vessel(3).backpropagate(omegas, ves6);

bi4_10_5 = bi4_10_5.backpropagate(omegas, ves5);

bi3_9_4 = bi3_9_4.backpropagate(omegas, bi4_10_5);

bi2_8_3 = bi2_8_3.backpropagate(omegas, bi3_9_4);

ves1 = ves1.backpropagate(omegas, bi2_8_3);

%% Forward Propagation

%Inlet of vessel 1
[Q1, P1] = ves1.forwardpropagate(omegas, ves1.s(0));

%Outlet of vessel 1
[Q1out, P1out] = ves1.forwardpropagate(omegas, ves1.s(ves1.L));

ves2 = bi2_8_3.vessel(1);

[Q2out, P2out] = ves2.forwardpropagate(omegas, ves2.s(ves2.L), P1out);

[Q8out, P8out, Q3out, P3out] = bi2_8_3.forwardpropagate(omegas, P2out);

[Q9out, P9out, Q4out, P4out] = bi3_9_4.forwardpropagate(omegas, P3out);

[Q10out, P10out, Q5out, P5out] = bi4_10_5.forwardpropagate(omegas, P4out);

[Q6out, P6out] = ves6.forwardpropagate(omegas, ves6.s(ves6.L), P5out);

[Q7out, P7out] = ves7.forwardpropagate(omegas, ves7.s(ves7.L), P6out);


%% Inverse Fourier Transform
sols = [Q1; P1; Q8out; P8out; Q9out; P9out; Q10out; P10out];

solt = InverseFourierTransform(t, omegas, sols, F);


[q1, p1, q11, p11, q12, p12] = deal(solt(1, :), solt(2, :), solt(3, :), solt(4, :), solt(5, :), solt(6, :));

[q13, p13] = deal(solt(7, :), solt(8, :));

xi = t;
%% Error calculations
%Average

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

emax13flow = max(abs((q13(2:end-1)-yilsuboutletflow(2:end-1)*10^-6)/max(yilsuboutletflow(2:end-1)*10^-6)));
emax13pressure = max(abs((p13-yilsuboutletpressure)./yilsuboutletpressure));

emax12flow = max(abs((q12(2:end-1)-yilccaoutletflow(2:end-1)*10^-6)/max(yilccaoutletflow(2:end-1)*10^-6)));
emax12pressure = max(abs((p12-yilccaoutletpressure)./yilccaoutletpressure));

emax11flow = max(abs((q11(2:end-1)-yibcoutletflow(2:end-1)*10^-6)/max(yibcoutletflow(2:end-1)*10^-6)));
emax11pressure = max(abs((p11-yibcoutletpressure)./yibcoutletpressure));

emax1pressure = max(abs((p1-yiinletpressure)./yiinletpressure));

%SYS

esys13flow = (max(q13(2:end-1))-max(yilsuboutletflow(2:end-1))*10^-6)/max(yilsuboutletflow(2:end-1)*10^-6);
esys13pressure = (max(p13)-max(yilsuboutletpressure))/max(yilsuboutletpressure);

esys12flow = (max(q12(2:end-1))-max(yilccaoutletflow(2:end-1))*10^-6)/max(yilccaoutletflow(2:end-1)*10^-6);
esys12pressure = (max(p12)-max(yilccaoutletpressure))/max(yilccaoutletpressure);

esys11flow = (max(q11(2:end-1))-max(yibcoutletflow(2:end-1)*10^-6))/max(yibcoutletflow(2:end-1)*10^-6);
esys11pressure = (max(p11)-max(yibcoutletpressure))/max(yibcoutletpressure);

esys1pressure = (max(p1)-max(yiinletpressure))/max(yiinletpressure);

%DIAS

edias13flow = (min(q13(2:end-1))-min(yilsuboutletflow(2:end-1))*10^-6)/max(yilsuboutletflow(2:end-1)*10^-6);
edias13pressure = (min(p13)-min(yilsuboutletpressure))/min(yilsuboutletpressure);

edias12flow = (min(q12(2:end-1))-min(yilccaoutletflow(2:end-1))*10^-6)/max(yilccaoutletflow(2:end-1)*10^-6);
edias12pressure = (min(p12)-min(yilccaoutletpressure))/min(yilccaoutletpressure);

edias11flow = (min(q11(2:end-1))-min(yibcoutletflow(2:end-1)*10^-6))/max(yibcoutletflow(2:end-1)*10^-6);
edias11pressure = (min(p11)-min(yibcoutletpressure))/min(yibcoutletpressure);

edias1pressure = (min(p1)-min(yiinletpressure))/min(yiinletpressure);

%% Plots
if (plotting)
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
end
