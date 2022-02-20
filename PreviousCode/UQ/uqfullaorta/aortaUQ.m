format long;
clear

% niPCE Parameters
ko = 2;    % chaos order
id = 10;    % # of uncertain variables
iPDF = [1 1 1 1 1 1 1 1 1 1]; %type of PDF (1)Gaussian (2)Uniform (3)Exponential
R1mean = 1; 
R1std = 0.05;
E1mean = 1; 
E1std = 0.05;
hmean = 1;
hstd = 0.05;
Lmean = 1;
Lstd = 0.05;
RW1mean = 1;
RW1std = 0.1;
RW2mean = 1;
RW2std = 0.1;
Cmean = 1;
Cstd = 0.1;
Qmean = 1;
Qstd = 0.1;
rhomean = 1;
rhostd = 0.02;
vmean = 1;
vstd = 0.02;

tic
[mean_L,std_L,l,t] = niPCE(ko,id,iPDF,E1mean,E1std,R1mean,R1std,hmean,hstd,Lmean,Lstd,RW1mean,RW1std,RW2mean,RW2std,Cmean,Cstd,Qmean,Qstd,rhomean,rhostd,vmean,vstd);
t1 = toc;

figure
fontxt = 21;
fontxt2 = 21;
k1 = [mean_L+std_L]*10^-3;
k2 = [mean_L-std_L]*10^-3;
% lightblue = [0.1,0.7,0.9];
lightgreen = [0.8,0.9,0.9];
patch([t' fliplr(t')], [k1 fliplr(k2)], lightgreen)
hold on
plot(t, (mean_L)*10^-3,'k','Linewidth',1.5)
%plot(t,(mean_L+std_L)*10^6,'b','linestyle','--','Linewidth',0.8)
%hold on
%plot(t,(mean_L-std_L)*10^6,'b','linestyle','--','Linewidth',0.8)
%hold on
box on
ax = gca;
ax.LineWidth = 2.5;
t1 = xlabel('time (s)');
t1.FontSize = fontxt;
t1 = ylabel('P (kPa)');
t1.FontSize = fontxt;
t1 = title('Uncertainty of P at r ren', 'FontSize', 21);
t1.FontSize = fontxt;
% t1 = legend('E+S','E');
% t1.FontSize = fontxt;
ax.XAxis.FontSize = fontxt;
ax.YAxis.FontSize = fontxt;


function [mean_Q,std_Q,l,t] = niPCE(ko,id,iPDF,E1mean,E1std,R1mean,R1std,hmean,hstd,Lmean,Lstd,RW1mean,RW1std,RW2mean,RW2std,Cmean,Cstd,Qmean,Qstd,rhomean,rhostd,vmean,vstd)
%=======Calculate Nodes and Weights from Tensor Product=====%
[sweights, snodes, psi, nhe, mxs] = GQprep1(ko,id,iPDF);
%====mean value and std of a=======
b=zeros(1+nhe,111);
%=======Integrate for spectral coefficients======
for ihq=1:mxs %mxs is number of evaluation nodes
    
    if mod(ihq,1000) == 0 ; disp(['niPCE ' num2str(ihq/mxs*100) ' % Completed']); end
    R1D = R1mean + R1std.*snodes(ihq,1); % Quadrature Node
    E1D = E1mean + E1std.*snodes(ihq,2);
    hD = hmean + hstd.*snodes(ihq,3);
    LD = Lmean + Lstd.*snodes(ihq,4);
    RW1D = RW1mean +RW1std.*snodes(ihq,5);
    RW2D = RW2mean +RW2std.*snodes(ihq,6);
    CD = Cmean + Cstd.*snodes(ihq,7);
    QD = Qmean + Qstd.*snodes(ihq,8);
    rhoD = rhomean + rhostd.*snodes(ihq,9);
    vD = vmean + vstd.*snodes(ihq,10);
    [Qoutlet,l,t] = aorta2(R1D,E1D,hD,LD,RW1D,RW2D,CD,QD,rhoD,vD);
    %Q(ihq) =  Qoutlet;
    for k = 1 : length(Qoutlet)
        
        for kk=1:nhe+1
            b(kk,k)=b(kk,k)+sweights(ihq)*psi(kk,ihq)*Qoutlet(k);
        end
        
    end
    
end

mean_Q(1:length(Qoutlet)) = b(1,1:length(Qoutlet));
std_Q(1:length(Qoutlet)) = zeros(1,length(Qoutlet));
 for ii=2:nhe+1
     std_Q(1:length(Qoutlet)) = std_Q(1:length(Qoutlet)) + b(ii,1:length(Qoutlet)).^2;
 end
 std_Q = std_Q.^0.5; 
end

