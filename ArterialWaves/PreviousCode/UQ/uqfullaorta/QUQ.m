format long;
clear

% niPCE Parameters
ko = 2;    % chaos order
id = 13;    % # of uncertain variables
iPDF = [1 1 1 1 1 1 1 1 1 1 1 1 1]; %type of PDF (1)Gaussian (2)Uniform (3)Exponential
ptmean = 1; 
ptstd = 0.1;
p1mean = 1; 
p1std = 0.1;
p2mean = 1;
p2std = 0.1;
p3mean = 1;
p3std = 0.1;
p4mean = 1;
p4std = 0.1;
p5mean = 1;
p5std = 0.1;
p6mean = 1;
p6std = 0.1;
p7mean = 1;
p7std = 0.1;
p8mean = 1;
p8std = 0.1;
p9mean = 1;
p9std = 0.1;
p10mean = 1;
p10std = 0.1;
p11mean = 1;
p11std = 0.1;
p12mean = 1;
p12std = 0.1;
outputs = 7;

load('../../automeris/aorta/BC.csv')

x = BC(1:end,1);
y = BC(1:end,2)*10^-6;
N = length(y);
T = BC(end,1);
xi = (x(1):T/((N)):T)';
yi = interp1q(x,y,xi);
Qin = yi(1:end);
t1 = xi(1:end);
N = length(Qin);

int = round(length(Qin)/12);
bounds = 13;
k(1) = 1;
for i = 2:bounds
    
    if i < bounds
        
        k(i) = k(i-1) + int;
        t{i-1} = t1(k(i-1):k(i),1);
        
    else
        
        k(i) = length(Qin);
        t{i-1} = t1(k(i-1):end,1);
        
    end
end

figure(1)
for i = 2:bounds
    
   
    plot(t{i-1},Qin(k(i-1):k(i))*10^6,'LineWidth',2)
    hold on
    
end

for i = 1:bounds-1
   
    l{i} = ['interval',' ',num2str(i),' : ',num2str(round(t{i}(1),4)),'-',num2str(round(t{i}(end),4)),' (s)'];
    
end

ld = 2.5;
colour1 = 'b';
colour2 = 'r';
fontxt = 14;
fontxt2 = 21;
t2 = legend(l);
t2.FontSize = fontxt;
t2 = xlabel('t(s)');
t2.FontSize = fontxt2;
t2 = ylabel('Q(ml/s)');
t2.FontSize = fontxt2;
box on
ax = gca;
ax.LineWidth = 2.5;
ax.XAxis.FontSize = fontxt2;
ax.YAxis.FontSize = fontxt2;


tic
[mean_L,std_L,l,t,b,nbox] = niPCE(ko,id,iPDF,ptmean,ptstd,p1mean,p1std,p2mean,p2std,p3mean,p3std,p4mean,p4std,p5mean,p5std,p6mean,p6std,p7mean,p7std,p8mean,p8std,p9mean,p9std,p10mean,p10std,p11mean,p11std,p12mean,p12std,outputs);
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

X = categorical({'Time-noise','interval 1','interval 2','interval 3','interval 4','interval 5','interval 6','interval 7','interval 8','interval 9','interval 10','interval 11','interval 12'});
X = reordercats(X,{'Time-noise','interval 1','interval 2','interval 3','interval 4','interval 5','interval 6','interval 7','interval 8','interval 9','interval 10','interval 11','interval 12'});

names = {'Brachiocephalic','L com. carotid','L subclavian','Sup.mesenteric','R renal','Inf mesenteric','R com. iliac'};
for i =1:outputs
    
    figure(i+1)
    bar(X,S(i,1:end))
    title('Sobol indices for Q at',names(i))
    xlabel('Uncertanties')
    ylabel('Percentage')

end



function [mean_Q,std_Q,l,t,b,nbox] = niPCE(ko,id,iPDF,ptmean,ptstd,p1mean,p1std,p2mean,p2std,p3mean,p3std,p4mean,p4std,p5mean,p5std,p6mean,p6std,p7mean,p7std,p8mean,p8std,p9mean,p9std,p10mean,p10std,p11mean,p11std,p12mean,p12std,outputs)
%=======Calculate Nodes and Weights from Tensor Product=====%
[sweights, snodes, psi, nhe, mxs,nbox] = GQprep1(ko,id,iPDF);
%====mean value and std of a=======
b=zeros(1+nhe,outputs);
%=======Integrate for spectral coefficients======
mxs
for ihq=1:mxs %mxs is number of evaluation nodes
    if mod(ihq,1000) == 0 ; disp(['niPCE ' num2str(ihq/mxs*100) ' % Completed']); end
    ptD = ptmean + ptstd.*snodes(ihq,1); % Quadrature Node
    p1D = p1mean + p1std.*snodes(ihq,2);
    p2D = p2mean + p2std.*snodes(ihq,3);
    p3D = p3mean + p3std.*snodes(ihq,4);
    p4D = p4mean + p4std.*snodes(ihq,5);
    p5D = p5mean +p5std.*snodes(ihq,6);
    p6D = p6mean +p6std.*snodes(ihq,7);
    p7D = p7mean + p7std.*snodes(ihq,8);
    p8D = p8mean + p8std.*snodes(ihq,9);
    p9D = p9mean + p9std.*snodes(ihq,10);
    p10D = p10mean + p10std.*snodes(ihq,11);
    p11D = p11mean + p11std.*snodes(ihq,12);
    p12D = p12mean + p12std.*snodes(ihq,13);
    [Qoutlet1,Qoutlet2,Qoutlet3,Qoutlet4,Qoutlet5,Qoutlet6,Qoutlet7,l,t] = aorta5(p1D,p2D,p3D,p4D,p5D,p6D,p7D,p8D,p9D,p10D,p11D,p12D,ptD);

    time1 = [1:length(Qoutlet1)];
    time2 = [1:length(Qoutlet1)];
    Q(1) =  max(Qoutlet1(time1));
    Q(2) =  max(Qoutlet2(time1));
    Q(3) =  max(Qoutlet3(time1));
    Q(4) =  max(Qoutlet4(time2));
    Q(5) =  max(Qoutlet5(time2));
    Q(6) =  max(Qoutlet6(time2));
    Q(7) =  max(Qoutlet7(time2));
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