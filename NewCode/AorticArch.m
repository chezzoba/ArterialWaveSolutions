%% Programme for implementation of the aortic arch model
clc;
clear;
addpath('../VesselModels/');

%% Determining Parameters

rho = 1050; % kg/m^3
PD = 1;

data = readtable('../Data/rsfs20170006supp1.xls', 'Sheet', sprintf('PD%d', PD));
As = table2array(readtable('../Data/rsfs20170006supp2.xls', 'Sheet', sprintf('PD%d', PD)));
Qs = table2array(readtable('../Data/rsfs20170006supp3.xls', 'Sheet', sprintf('PD%d', PD)));
[As, Qs] = deal(As(:, 2:end), Qs(:, 2:end));

Nx = height(Qs);

avQ = mean(Qs, 2);
Qdrops = find(abs((avQ(1:Nx-1) - avQ(2:Nx)) ./ avQ(1:Nx-1)) > 0.1);
pos1 = Qdrops(1);
pos3 = Qdrops(end);

poses = round([1, pos1/2, pos1, (pos1 + pos3) / 2,...
    pos3, (Nx + 2*pos3)/3, (2*Nx + pos3)/3, Nx]);

xs = data.AnalysisPlane(poses) * 1e-3;
Asysp = data.MaximumAreaMeasured(poses) * 1e-6;
Adiasp = data.MinimumAreaMeasured(poses) * 1e-6;
Psysp = data.MaximumPressure(poses) * 133.32237;
Pdiasp = data.MinimumPressure(poses) * 133.32237;

Ls = (xs(2:end) - xs(1:length(xs)-1));
Rs = data.EndDiastolicRadius(poses) * 1e-3;
betas = (1 - sqrt(Adiasp ./ Asysp)) ./ (Rs .* (Psysp - Pdiasp));
as = atan((Rs(1:length(Rs)-1) - Rs(2:end)) ./ Ls);

parents = [2, 3, 4];
chila = [3, 4, 5];
Rpparents = sqrt(2 .* rho ./ betas(parents)) ./ (pi .* Rs(parents) .^ 2.5);
Rpchila = sqrt(2 .* rho ./ betas(chila)) ./ (pi .* Rs(chila) .^ 2.5);
Rpsa = (Rpparents .* Rpchila) ./ (Rpparents - Rpchila);
betasa = 2 * rho ./ (Rpsa .* pi .* Rs(chila) .^ 2.5) .^ 2;
Lsa = 2e-2;

plotting = true;

%% Importing the data from the Flores plots

WKP7 = [17513091.0464560, 214342021.080118, 4.67859041979278e-09];

ves1 = Vessel(Rs(1), Ls(1), as(1), betas(1),...
    rho, [0, 0, 0]);
ves1.type = 5;

bi2_8_3 = Bifurcation([Rs(2), Rs(3), Rs(3)], [Ls(2), Lsa, Ls(3)],...
    [betas(2), betasa(1), betas(3)],...
    rho, [as(1), 0.001, as(3)],...
    51918000, 1060800000, 8.6974e-10);

bi3_9_4 = Bifurcation([Rs(3) Rs(4) Rs(4)], [Ls(3) Lsa Ls(4)],...
    [betas(3),betasa(2),betas(4)], rho,...
    [as(3),1e-03,as(4)],...
    191515000, 5.22129e+09, 1.767e-10);

bi4_10_5 = Bifurcation([Rs(4),Rs(5),Rs(5)], [Ls(4),Lsa,Ls(5)],...
    [betas(4),betasa(3),betas(5)], rho,...
    [as(4),1e-03,as(5)],...
    98820000, 1.30183e+09, 7.0871e-10);

ves6 = Vessel(Rs(6), Ls(6), as(6), betas(6), rho, [0, 0, 0]);

ves7 = Vessel(Rs(7), Ls(7), as(7), betas(7), rho, WKP7);
ves7.type = 2;

BC = [(0:24) / 25; Qs(1, :)].';
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
sols = [Q1; P1; Q8out; P8out; Q9out; P9out; Q10out; P10out; Q7out; P7out];

solt = InverseFourierTransform(t, omegas, sols, F);


[q1, p1, q8, p8, q9, p9] = deal(solt(1, :), solt(2, :), solt(3, :), solt(4, :), solt(5, :), solt(6, :));

[q10, p10, q7, p7] = deal(solt(7, :), solt(8, :), solt(9, :), solt(10, :));



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
    
    figure
    plot(t,p1*10^-3,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2 = title('Pressure - inlet (1)');
    t2.FontSize = fontxt2;
    t2 = xlabel('time (s)');
    t2.FontSize = fontxt2;
    t2 = ylabel('P (kPa)');
    t2.FontSize = fontxt2;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.Layer = 'top';
    ax.XAxis.FontSize = fontxt2;
    ax.YAxis.FontSize = fontxt2;
    
    figure
    plot(t,q8*10^6,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2 = title('Volume flow rate - bc outlet (8)');
    t2.FontSize = fontxt2;
    t2 = xlabel('time (s)');
    t2.FontSize = fontxt2;
    t2 = ylabel('Q (ml/s)');
    t2.FontSize = 13;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.Layer = 'top';
    ax.XAxis.FontSize = fontxt2;
    ax.YAxis.FontSize = fontxt2;
    
    figure
    plot(t,p8*10^-3,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2.FontSize = fontxt2;
    t2 = title('Pressure -  bc outlet (8)');
    t2.FontSize = fontxt2;
    t2 = xlabel('time (s)');
    t2.FontSize = fontxt2;
    t2 = ylabel('P (kPa)');
    t2.FontSize = fontxt2;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.Layer = 'top';
    ax.XAxis.FontSize = fontxt2;
    ax.YAxis.FontSize = fontxt2;
    
    figure
    plot(t,q9*10^6,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2.FontSize = fontxt2;
    t2 = title('Volume flow rate - lcca outlet (9)');
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
    
    figure
    plot(t,p9*10^-3,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2.FontSize = fontxt2;
    t2 = title('Pressure - lcca outlet (9)');
    t2.FontSize = fontxt2;
    t2 = xlabel('time (s)');
    t2.FontSize = fontxt2;
    t2 = ylabel('P (kPa)');
    t2.FontSize = fontxt2;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.Layer = 'top';
    ax.XAxis.FontSize = fontxt2;
    ax.YAxis.FontSize = fontxt2;
    
    figure
    plot(t,q10*10^6,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2.FontSize = fontxt2;
    t2 = title('Volume flow rate - lsub outlet (10)');
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

    figure
    plot(t,p10*10^-3,colour1,'LineWidth',ld)
    hold on
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2.FontSize = fontxt2;
    t2 = title('Pressure - lcca outlet (10)');
    t2.FontSize = fontxt2;
    t2 = xlabel('time (s)');
    t2.FontSize = fontxt2;
    t2 = ylabel('P (kPa)');
    t2.FontSize = fontxt2;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.Layer = 'top';
    ax.XAxis.FontSize = fontxt2;
    ax.YAxis.FontSize = fontxt2;
    
    figure
    plot(t,q7*10^6,colour1,'LineWidth',ld)
    hold on
    plot((0:24) / 25, Qs(end, :),colour2,'LineWidth',ld)
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2.5;
    t2 = legend('1-D (Present)','4-D');
    t2.FontSize = fontxt2;
    t2 = title('Volume flow rate - outlet (7)');
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
end
