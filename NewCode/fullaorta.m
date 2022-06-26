%% Programme for implementation of the full aorta model
clc;
clear;
addpath('../VesselModels/');

plotting = false;

tic
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
be = (1 - v^2) ./ (E .* h);

%% Importing the data from the Flores plots
load('../PreviousCode/automeris/aorta/BC.csv')
load('../PreviousCode/automeris/aorta/Inlet_Pressure.csv')
load('../PreviousCode/automeris/aorta/bc_outlet_11_flow.csv')
load('../PreviousCode/automeris/aorta/bc_outlet_11_Pressure.csv')
load('../PreviousCode/automeris/aorta/lcca_outlet_12_flow.csv')
load('../PreviousCode/automeris/aorta/lcca_outlet_12_Pressure.csv')
load('../PreviousCode/automeris/aorta/lsub_outlet_13_flow.csv')
load('../PreviousCode/automeris/aorta/lsub_outlet_13_Pressure.csv')
load('../PreviousCode/automeris/aorta/sma_outlet_15_flow.csv')
load('../PreviousCode/automeris/aorta/sma_outlet_15_Pressure.csv')
load('../PreviousCode/automeris/aorta/rena_outlet_16_flow.csv')
load('../PreviousCode/automeris/aorta/rena_outlet_16_Pressure.csv')
load('../PreviousCode/automeris/aorta/imma_outlet_18_flow.csv')
load('../PreviousCode/automeris/aorta/imma_outlet_18_Pressure.csv')
load('../PreviousCode/automeris/aorta/riliac_outlet_19_flow.csv')
load('../PreviousCode/automeris/aorta/riliac_outlet_19_Pressure.csv')



bi10_19_20 = Bifurcation([R(10), R(19), R(20)], [L(10), L(19), L(20)],...
    [be(10), be(19), be(20)], rho, [a(10), a(19), a(20)],...
    [RW1(19), RW1(20)], [RW2(19), RW2(20)], [CWK(19), CWK(20)]);

bi10_19_20.type = 2;

bi9_18_10 = Bifurcation([R(9), R(18), R(10)], [L(9), L(18), L(10)],...
    [be(9), be(18), be(10)], rho, [a(9), a(18), a(10)],...
    RW1(18), RW2(18), CWK(18));

bi8_17_9 = Bifurcation([R(8), R(17), R(9)], [L(8), L(17), L(9)],...
    [be(8), be(17), be(9)], rho, [a(8), a(17), a(9)],...
    RW1(17), RW2(17), CWK(17));

bi7_16_8 = Bifurcation([R(7), R(16), R(8)], [L(7), L(16), L(8)],...
    [be(7), be(16), be(8)], rho, [a(7), a(16), a(8)],...
    RW1(16), RW2(16), CWK(16));

bi6_15_7 = Bifurcation([R(6), R(15), R(7)], [L(6), L(15), L(7)],...
    [be(6), be(15), be(7)], rho, [a(6), a(15), a(7)],...
    RW1(15), RW2(15), CWK(15));

bi5_14_6 = Bifurcation([R(5), R(14), R(6)], [L(5), L(14), L(6)],...
    [be(5), be(14), be(6)], rho, [a(5), a(14), a(6)],...
    RW1(14), RW2(14), CWK(14));

bi3_13_4 = Bifurcation([R(3), R(13), R(4)], [L(3), L(13), L(4)],...
    [be(3), be(13), be(4)], rho, [a(3), a(13), a(4)],...
    RW1(13), RW2(13), CWK(13));

ves4 = bi3_13_4.vessel(3);

bi2_12_3 = Bifurcation([R(2), R(12), R(3)], [L(2), L(12), L(3)],...
    [be(2), be(12), be(3)], rho, [a(2), a(12), a(3)],...
    RW1(12), RW2(12), CWK(12));

bi1_11_2 = Bifurcation([R(1), R(11), R(2)], [L(1), L(11), L(2)],...
    [be(1), be(11), be(2)], rho, [a(1), a(11), a(2)],...
    RW1(11), RW2(11), CWK(11));
bi1_11_2.type = 5;

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


%% Calculations in the freequency domain
omegas = om * (0:nh);
omegas(1) = 1e-10;


   
%% Backward Propagation  

%10-19-20 --> Type II Bifurcation
bi10_19_20 = bi10_19_20.backpropagate(omegas);

%9-18-10 --> Type III Bifurcation
bi9_18_10 = bi9_18_10.backpropagate(bi10_19_20);

%8-17-9 --> Type III Bifurcation
bi8_17_9 = bi8_17_9.backpropagate(bi9_18_10);

%7-16-8 --> Type III Bifurcation
bi7_16_8 = bi7_16_8.backpropagate(bi8_17_9);

%6-15-7 --> Type III Bifurcationzi
bi6_15_7 = bi6_15_7.backpropagate(bi7_16_8);

%5-14-6 --> Type III Bifurcation
bi5_14_6 = bi5_14_6.backpropagate(bi6_15_7);

%4-5 --> Two connected vessels
ves4 = ves4.backpropagate(bi5_14_6);

%3-13-4 --> Type III Bifurcation  
bi3_13_4 = bi3_13_4.backpropagate(ves4);

%2-12-3 --> Type III Bifurcation    
bi2_12_3 = bi2_12_3.backpropagate(bi3_13_4);

%1-11-2 --> Type V Bifurcation 
bi1_11_2 = bi1_11_2.backpropagate(bi2_12_3);

%% Forward Propagation 
%Inlet of vessel 1
ves1 = bi1_11_2.vessel(1);

[Q1, P1] = ves1.forwardpropagate(ves1.s(0), 0);

%Outlet of vessel 1
[Q1out, P1out] = ves1.forwardpropagate(ves1.s(ves1.L), 0);

% Brachiocephalic and AO II
[Q11out, P11out, Q2out, P2out] = bi1_11_2.forwardpropagate(P1out);

%L com. carotid and AO III
[Q12out, P12out, Q3out, P3out] = bi2_12_3.forwardpropagate(P2out);

%Left subclavian and AO IV
[Q13out, P13out, Q4out, P4out] = bi3_13_4.forwardpropagate(P3out);

%AO V
ves5 = bi5_14_6.vessel(1);
[Q5out, P5out] = ves5.forwardpropagate(ves5.s(ves5.L), P4out);

%AO VI
[Q14out, P14out, Q6out, P6out] = bi5_14_6.forwardpropagate(P5out);

%Sup. mesenteric and AO VII
[Q15out, P15out, Q7out, P7out] = bi6_15_7.forwardpropagate(P6out);

%Renal and AO VIII
[Q16out, P16out, Q8out, P8out] = bi7_16_8.forwardpropagate(P7out);

%AO IX and 17
[Q17out, P17out, Q9out, P9out] = bi8_17_9.forwardpropagate(P8out);

%Inf. mesenteric and AO X
[Q18out, P18out, Q10out, P10out] = bi9_18_10.forwardpropagate(P9out);

%R com. iliac
[Q19out, P19out, Q20out, P20out] = bi10_19_20.forwardpropagate(P10out);


% Inverse Fourier Transform
sols = [Q1; P1; Q11out; P11out; Q12out; P12out; Q13out; P13out;
    Q15out; P15out; Q16out; P16out; Q18out; P18out; Q19out; P19out];

solt = InverseFourierTransform(t, -omegas, conj(sols), F);

[q1, p1, q11, p11, q12, p12] = deal(solt(1, :), solt(2, :), solt(3, :), solt(4, :), solt(5, :), solt(6, :));
[q13, p13, q15, p15, q16, p16] = deal(solt(7, :), solt(8, :), solt(9, :), solt(10, :), solt(11, :), solt(12, :));
[q18, p18, q19, p19] = deal(solt(13, :), solt(14, :), solt(15, :), solt(16, :));

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
end
