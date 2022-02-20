format long;
clear
tic
% niPCE Parameters
ko = 1;    % chaos order
id = 20;    % # of uncertain variables
iPDF = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %type of PDF (1)Gaussian (2)Uniform (3)Exponential
Rinmean = [15.2,13.9,13.7,13.5,12.3,9.9,9.7,9.62,9.55,9.07,6.35,3.6,4.8,4.45,3.75,2.8,2.8,2,6,6]*10^-3;
Rinstd = Rinmean*0.1;
outputs = 1;    %Number of code outputs (QoI)

tic
[mean_L,std_L,l,t,b,nbox] = niPCE(ko,id,iPDF,outputs,Rinmean,Rinstd);
toc

%Create a matrix with information about the uncertainties that participated
%in the calculation of each coefficient. The nbox from the GQprep function
%is required.
for ind = 1:id 
    r = 0;

    for i = 1:length(nbox(1:end,1))

         k = 0;

         if nbox(i,ind) ~= 0

             k = k + 1;

         end

         if k ~= 0


            r = r + 1;
            row(ind,r) = i;

        else

        end

    end

end

%Calculate the participation of each uncertainty to the total std
S = zeros(outputs,id);
for i = 1:id

    for j = 1:length(row(1,1:end))
       
        for k = 1:outputs
            
            S(k,i) = S(k,i) + b(row(i,j),k)^2;
        
        end
    end

end

s = std_L.^2;

%Divide S with the standard deviation in order to calculate the
%index in precentage
for i = 1:outputs
    
    S(i,1:end) = S(i,1:end)./s(i);
    
end

X = categorical({'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R17','R18','R19','R20'});
X = reordercats(X,{'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16','R17','R18','R19','R20'});

names = {'Brachiocephalic (11)','L com. carotid (12)','L subclavian (13)','Sup.mesenteric (15)','R renal (16)','Inf mesenteric (18)','R com. iliac (19)'};
for i =1:outputs
    
    figure(i)
    bar(X,S(i,1:end))
    title('Sobol indices for systolic Q at',names(i))
    xlabel('Uncertanties')
    ylabel('Percentage')

end

toc


function [mean_Q,std_Q,l,t,b,nbox] = niPCE(ko,id,iPDF,outputs,Rinmean,Rinstd)
%=======Calculate Nodes and Weights from Tensor Product=====%
[sweights, snodes, psi, nhe, mxs,nbox] = GQprep1(ko,id,iPDF);
%====mean value and std of a=======
b=zeros(1+nhe,outputs);
mxs   %Print the number of code evaluations
%=======Integrate for spectral coefficients======
for ihq=1:mxs %mxs is number of evaluation nodes
    if mod(ihq,5000) == 0 ; disp(['niPCE ' num2str(ihq/mxs*100) ' % Completed']); end
    for i = 1 :id
       
        RinD(i) = Rinmean(i) + Rinstd(i).*snodes(ihq,i); % Quadrature Node
        
    end
    [Qoutlet1,l,t] = aorta3(RinD);
%     
%   time1 = [3:40];                                 %Systolic Phase
%   time2 = [10:46];                                %Systolic Phase
%   time1 = [1:2,41:length(Qoutlet1)];              %Diastolic Phase
%   time2 = [1:9,47:length(Qoutlet1)];              %Diastolic Phase
  time1 = [1:length(Qoutlet1)];                     %Full Cardiac Cycle
  time2 = [1:length(Qoutlet1)];                     %Full Cardiac Cycle

    Q(1) =  Qoutlet1(time1);
%     Q(2) =  Qoutlet2(time1);
%     Q(3) =  Qoutlet3(time1);
%     Q(4) =  Qoutlet4(time2);
%     Q(5) =  Qoutlet5(time2);
%     Q(6) =  Qoutlet6(time2);
%     Q(7) =  Qoutlet7(time2);
    for k = 1:outputs
        
        for kk=1:nhe+1
            b(kk,k)=b(kk,k)+sweights(ihq)*psi(kk,ihq)*Q(k);
        end
        
    end
    
end

mean_Q(1,1:outputs) = b(1,1:outputs);
std_Q = zeros(1,outputs);
 for ii=2:nhe+1
     std_Q = std_Q +b(ii,1:outputs).^2;
 end
 std_Q = std_Q.^0.5; 
end