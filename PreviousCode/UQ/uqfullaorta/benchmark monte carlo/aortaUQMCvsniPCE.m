format long;
clear

% niPCE Parameters
% niPCE Parameters
ko = 1;    % chaos order
id = 20;    % # of uncertain variables
iPDF = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %type of PDF (1)Gaussian (2)Uniform (3)Exponential
Rinmean = [15.2,13.9,13.7,13.5,12.3,9.9,9.7,9.62,9.55,9.07,6.35,3.6,4.8,4.45,3.75,2.8,2.8,2,6,6]*10^-3;
Rinstd = Rinmean*0.1;
Samples = 5000;
outputs = 7;

tic
[mL,sL,Jmean,Jstd,J,l,t] = MonteCarlo(Rinmean,Rinstd,Samples);
t2 = toc;

tic
[mean_L,std_L,l,t] = niPCE(ko,id,iPDF,Rinmean,Rinstd);
t1 = toc;

figure
plot(t, mean_L)
hold on
plot(t,mean_L+std_L)
hold on
plot(t,mean_L-std_L)
hold on
xlabel('time (s)')
ylabel('Q (ml/s)')
hold on
plot(t, mL,'--')
hold on
plot(t,mL+sL,'--')
hold on
plot(t,mL-sL,'--')
xlabel('time(s)')
ylabel('P(ml/s)')
title('Uncertainty of Flow at the outlet MC(N=50000) vs niPCE(C=1)')



function [mean_Q,std_Q,l,t] = niPCE(ko,id,iPDF,Rinmean,Rinstd)
%=======Calculate Nodes and Weights from Tensor Product=====%
[sweights, snodes, psi, nhe, mxs] = GQprep1(ko,id,iPDF);
%====mean value and std of a=======
b=zeros(1+nhe,111);
%=======Integrate for spectral coefficients======
mxs
for ihq=1:mxs %mxs is number of evaluation nodes
 
    for i = 1 :id
       
        RinD(i) = Rinmean(i) + Rinstd(i).*snodes(ihq,i); % Quadrature Node
        
    end
    
    [Qoutlet,l,t] = aorta3(RinD);
 
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



function [mL,sL,Jmean,Jstd,J,l,t] = MonteCarlo(Rinmean,Rinstd,N)

tic
%J = zeros(N,1); Jmean = J; Jstd = J; % Preallocation
for ii=1:N
    
    for j = 1 : length(Rinmean)
        
        RinD(j) = normrnd(Rinmean(j),Rinstd(j));
        
    end
    
    [J1,l,t] = aorta3(RinD);
    J(ii,1:length(J1)) = J1;
    for k = 1 : length(J1)
        
        Jmean(ii,k) = mean(J(1:ii,k)); Jstd(ii,k) = std(J(1:ii,k));
        
    end   
        if mod(ii,500) == 0 ; disp(['Monte-Carlo ' num2str(ii/N*100) ' % Completed']); end
end           
toc
mL(1:length(J1)) = mean(J(1:end,1:k));
sL(1:length(J1)) =  std(J(1:end,1:k));
end

