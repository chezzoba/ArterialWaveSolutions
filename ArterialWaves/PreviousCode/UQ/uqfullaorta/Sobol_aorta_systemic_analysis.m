format long;
clear

prompt = 'Please enter: 1 for the full cardiac cycle, 2 for the systolic phase, 3 for the diastolic phase, 4 for the maximum value ';
x = input(prompt);

% niPCE Parameters
ko = 2;    % chaos order
id = 10;    % # of uncertain variables
iPDF = [1 1 1 1 1 1 1 1 1 1]; %type of PDF (1)Gaussian (2)Uniform (3)Exponential
Rmean = 1; 
Rstd = 0.1;
Emean = 1; 
Estd = 0.1;
hmean = 1;
hstd = 0.1;
vmean = 1;
vstd = 0.1;
Lmean = 1;
Lstd = 0.1;
RW1mean = 1;
RW1std = 0.1;
RW2mean = 1;
RW2std = 0.1;
Cmean = 1;
Cstd = 0.1;
Qmean = 1;
Qstd = 0.1;
rhomean = 1;
rhostd = 0.1;
outputs = 7;

tic
[mean_L,std_L,l,t,b,nbox] = niPCE(ko,id,iPDF,Lmean,Lstd,RW1mean,RW1std,RW2mean,RW2std,Cmean,Cstd,Emean,Estd,Rmean,Rstd,hmean,hstd,vmean,vstd,rhomean,rhostd,Qmean,Qstd,outputs,x);
toc

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

S = zeros(outputs,id);
for i = 1:id

    for j = 1:length(row(1,1:end))
       
        for k = 1:outputs
            
            S(k,i) = S(k,i) + b(row(i,j),k)^2;
        
        end
    end

end

s = std_L.^2;

for i = 1:outputs
    
    S(i,1:end) = S(i,1:end)./s(i);
    
end

X = categorical({'Q','R','E','L','H','rho' , 'v','RW1','RW2','C'});
X = reordercats(X,{'Q','R','E','L','H','rho' , 'v','RW1','RW2','C'});

names = {'Brachiocephalic','L com. carotid','L subclavian','Sup.mesenteric','R renal','Inf mesenteric','R com. iliac'};
for i =1:outputs
    
    figure(i)
    bar(X,S(i,1:end))
    title('Sobol indices for Q at',names(i))
    xlabel('Uncertanties')
    ylabel('Percentage')

end


function [mean_Q,std_Q,l,t,b,nbox] = niPCE(ko,id,iPDF,Lmean,Lstd,RW1mean,RW1std,RW2mean,RW2std,Cmean,Cstd,Emean,Estd,Rmean,Rstd,hmean,hstd,vmean,vstd,rhomean,rhostd,Qmean,Qstd,outputs,x)
%=======Calculate Nodes and Weights from Tensor Product=====%
[sweights, snodes, psi, nhe, mxs,nbox] = GQprep1(ko,id,iPDF);
%====mean value and std of a=======
b=zeros(1+nhe,outputs);
mxs
%=======Integrate for spectral coefficients======
for ihq=1:mxs %mxs is number of evaluation nodes
    if mod(ihq,500) == 0 ; disp(['niPCE ' num2str(ihq/mxs*100) ' % Completed']); end
    QD = Qmean + Qstd.*snodes(ihq,1);
    RD = Rmean + Rstd.*snodes(ihq,2); % Quadrature Node
    ED = Emean + Estd.*snodes(ihq,3);
    LD = Lmean + Lstd.*snodes(ihq,4);
    hD = hmean + hstd.*snodes(ihq,5);
    rhoD = rhomean + rhostd.*snodes(ihq,6);
    vD = vmean + vstd.*snodes(ihq,7);
    RW1D = RW1mean +RW1std.*snodes(ihq,8);
    RW2D = RW2mean +RW2std.*snodes(ihq,9);
    CD = Cmean + Cstd.*snodes(ihq,10);
    
    [Qoutlet1,Qoutlet2,Qoutlet3,Qoutlet4,Qoutlet5,Qoutlet6,Qoutlet7,l,t] = aorta(RD,ED,hD,LD,RW1D,RW2D,CD,QD,rhoD,vD);
        
%     Uncomment to select the area of interest
    
    if x == 1 || x == 4
        time1 = [1:length(Qoutlet1)];               %mean
        time2 = [1:length(Qoutlet1)];               %mean
    elseif x == 2
        time1 = [3:40];                           %systolic
        time2 = [10:46];                          %systolic
    elseif x ==3 
        time1 = [1:2,41:length(Qoutlet1)];        %diastolic
        time2 = [1:9,47:length(Qoutlet1)];        %diastolic
    else
        
        disp('Worng input! Please sepcify one of the available options')
        return
    end
    
    if x == 1 || x == 2 || x == 3
        Q(1) =  mean(Qoutlet1(time1));
        Q(2) =  mean(Qoutlet2(time1));
        Q(3) =  mean(Qoutlet3(time1));
        Q(4) =  mean(Qoutlet4(time2));
        Q(5) =  mean(Qoutlet5(time2));
        Q(6) =  mean(Qoutlet6(time2));
        Q(7) =  mean(Qoutlet7(time2));
        
    elseif x == 4
        Q(1) =  max(Qoutlet1(time1));
        Q(2) =  max(Qoutlet2(time1));
        Q(3) =  max(Qoutlet3(time1));
        Q(4) =  max(Qoutlet4(time2));
        Q(5) =  max(Qoutlet5(time2));
        Q(6) =  max(Qoutlet6(time2));
        Q(7) =  max(Qoutlet7(time2));
    end
    
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