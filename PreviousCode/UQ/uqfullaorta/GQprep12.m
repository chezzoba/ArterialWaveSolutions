%Kyriakos D. Kantarakias 15/11/18
%
%This script creates the data structure required for the PCE
%IMPORTANT!!! psi matrix is (1:nhe+1,mxs) and not (0:nhe,mxs)
%input:
%    id       = # of uncertain variables
%           ko       = chaos order
%           iPDF(id) = indicates the PDF of ith uncertain variable
%                      iPDF(i)=1 for normal distribution (Hermite)
%                      iPDF(i)=2 for uniform distribution (Legendre)
%                      iPDF(i)=3 for exponential distribution (Laguerre)
%Output:
%              nhe +1 = # of polynomials
%              mxs    = # of Smolyak nodes
%              psi(0:nhe,1:mxs) = ith polynomial at jth Gauss or Smolyak node
%              sweights(mxs) = Gauss or Smolyak quadrature weights
%
%              Output is saved in file 'GQ.dat'

function [sweights, snodes, psi, nhe, mxs,nbox] = GQprep1(ko,id,iPDF)
mxhq=20;


if id == 0
   disp('No uncertain variables are defined.')
   nhe=0;
   mxs=1;
   sweights(1)=1;
   psi(0,1)   =1;
   return
end
    nhe = nchoosek(ko+id,ko);
    ks  = ko+id;
    if ks > mxhq 
        disp('Increase order of univariate rules ks')
        return
    end
    mxind = nchoosek(ks-1,id-1);
    %ind   = zeros(mxind,id); do I need preallocation here?
    l1    = max(id,ks-id+1);
    mxs   = 0;
    for l= l1:ks
        [kount,ind] = drop(id,l,mxind);
        for i=1:kount
            ldim = 1; 
            for j=1:id
                ldim=ldim*ind(i,j);
            end    
            mxs = mxs + ldim;
        end
    end 
    nhe = nhe - 1;    
    
    uweights = zeros(mxhq,mxhq,id);
    unodes = zeros(mxhq,mxhq,id);

%     =========
      for i=1:id
%     =========
         if iPDF(i)== 1 %Hermite
%
            unodes  (1,1,i)=   0; 
            uweights(1,1,i)=   1.772453850906;
%
            unodes  (2,1,i)=   0.707106781186548;
            unodes  (2,2,i)= - 0.707106781186548;
            uweights(2,1,i)=   0.8862269254528;
            uweights(2,2,i)=   0.8862269254528;
%
            unodes  (3,1,i)=   0;
            unodes  (3,2,i)=   1.224744871391589;
            unodes  (3,3,i)= - 1.224744871391589;
            uweights(3,1,i)=   1.181635900604;
            uweights(3,2,i)=   0.2954089751509;
            uweights(3,3,i)=   0.2954089751509;
%
            unodes  (4,1,i)=   0.524647623275290;
            unodes  (4,2,i)=   1.650680123885785;
            unodes  (4,3,i)= - 0.524647623275290;
            unodes  (4,4,i)= - 1.650680123885785;
            uweights(4,1,i)=   0.8049140900055;
            uweights(4,2,i)=   0.08131283544725;
            uweights(4,3,i)=   0.8049140900055;
            uweights(4,4,i)=   0.08131283544725;
%
            unodes  (5,1,i)=   0;
            unodes  (5,2,i)=   0.958572464613819;
            unodes  (5,3,i)=   2.020182870456086;
            unodes  (5,4,i)= - 0.958572464613819;
            unodes  (5,5,i)= - 2.020182870456086;
            uweights(5,1,i)=   0.9453087204829;
            uweights(5,2,i)=   0.3936193231522;
            uweights(5,3,i)=   0.01995324205905;
            uweights(5,4,i)=   0.3936193231522;
            uweights(5,5,i)=   0.01995324205905;
%
            unodes  (6,1,i)=   0.436077411927617;
            unodes  (6,2,i)=   1.335849074013697;
            unodes  (6,3,i)=   2.350604973674492;
            unodes  (6,4,i)= - 0.436077411927617;
            unodes  (6,5,i)= - 1.335849074013697;
            unodes  (6,6,i)= - 2.350604973674492;
            uweights(6,1,i)=   0.7246295952244;
            uweights(6,2,i)=   0.1570673203229;
            uweights(6,3,i)=   0.004530009905509;
            uweights(6,4,i)=   0.7246295952244;
            uweights(6,5,i)=   0.1570673203229;
            uweights(6,6,i)=   0.004530009905509;
%
            unodes  (7,1,i)=   0;
            unodes  (7,2,i)=   0.816287882858965;
            unodes  (7,3,i)=   1.673551628767471;
            unodes  (7,4,i)=   2.651961356835233;
            unodes  (7,5,i)= - 0.816287882858965;
            unodes  (7,6,i)= - 1.673551628767471;
            unodes  (7,7,i)= - 2.651961356835233;
            uweights(7,1,i)=   0.8102646175568;
            uweights(7,2,i)=   0.4256072526101;
            uweights(7,3,i)=   0.05451558281913;
            uweights(7,4,i)=   0.0009717812450995;
            uweights(7,5,i)=   0.4256072526101;
            uweights(7,6,i)=   0.05451558281913;
            uweights(7,7,i)=   0.0009717812450995;
%
            unodes  (8,1,i)=   0.381186990207322;
            unodes  (8,2,i)=   1.157193712446780;
            unodes  (8,3,i)=   1.981656756695843;
            unodes  (8,4,i)=   2.930637420257244;
            unodes  (8,5,i)= - 0.381186990207322;
            unodes  (8,6,i)= - 1.157193712446780;
            unodes  (8,7,i)= - 1.981656756695843;
            unodes  (8,8,i)= - 2.930637420257244;
            uweights(8,1,i)=   0.6611470125582;
            uweights(8,2,i)=   0.2078023258149;
            uweights(8,3,i)=   0.01707798300741;
            uweights(8,4,i)=   0.0001996040722114;
            uweights(8,5,i)=   0.6611470125582;
            uweights(8,6,i)=   0.2078023258149;
            uweights(8,7,i)=   0.01707798300741;
            uweights(8,8,i)=   0.0001996040722114;
%
            unodes  (9,1,i)=   0;
            unodes  (9,2,i)=   0.723551018752838;
            unodes  (9,3,i)=   1.468553289216668;
            unodes  (9,4,i)=   2.266580584531843;
            unodes  (9,5,i)=   3.190993201781528;
            unodes  (9,6,i)= - 0.723551018752838;
            unodes  (9,7,i)= - 1.468553289216668;
            unodes  (9,8,i)= - 2.266580584531843;
            unodes  (9,9,i)= - 3.190993201781528;
            uweights(9,1,i)=   0.7202352156061;
            uweights(9,2,i)=   0.4326515590026;
            uweights(9,3,i)=   0.08847452739438;
            uweights(9,4,i)=   0.004943624275537;
            uweights(9,5,i)=   0.00003960697726326;
            uweights(9,6,i)=   0.4326515590026;
            uweights(9,7,i)=   0.08847452739438;
            uweights(9,8,i)=   0.004943624275537;
            uweights(9,9,i)=   0.00003960697726326;
%
            unodes  (10,1,i)=   0.342901327;
            unodes  (10,2,i)=   1.036610830;
            unodes  (10,3,i)=   1.756683649;
            unodes  (10,4,i)=   2.532731674;
            unodes  (10,5,i)=   3.436159119;
            unodes  (10,6,i)= - 0.342901327;
            unodes  (10,7,i)= - 1.036610830;
            unodes  (10,8,i)= - 1.756683649;
            unodes  (10,9,i)= - 2.532731674;
            unodes  (10,10,i)=- 3.436159119;
            uweights(10,1,i)=   0.610862634;
            uweights(10,2,i)=   0.240138611;
            uweights(10,3,i)=   0.0338743945;
            uweights(10,4,i)=   0.00134364577;
            uweights(10,5,i)=   0.00000764043286;
            uweights(10,6,i)=   0.610862634;
            uweights(10,7,i)=   0.240138611;
            uweights(10,8,i)=   0.0338743945;
            uweights(10,9,i)=   0.00134364577;
            uweights(10,10,i)=  0.00000764043286;
%
            unodes  (11,1 ,i)=   0;
            unodes  (11,2 ,i)=   0.656809566882100;
            unodes  (11,3 ,i)=   1.326557084494933;
            unodes  (11,4 ,i)=   2.025948015825755;
            unodes  (11,5 ,i)=   2.783290099781652;
            unodes  (11,6 ,i)=   3.668470846559583;
            unodes  (11,7 ,i)= - 0.656809566882100;
            unodes  (11,8 ,i)= - 1.326557084494933;
            unodes  (11,9 ,i)= - 2.025948015825755;
            unodes  (11,10,i)= - 2.783290099781652;
            unodes  (11,11,i)= - 3.668470846559583;
            uweights(11,1 ,i)=   0.6547592869146;
            uweights(11,2 ,i)=   0.4293597523561;
            uweights(11,3 ,i)=   0.1172278751677;
            uweights(11,4 ,i)=   0.01191139544491;
            uweights(11,5 ,i)=   0.0003468194663233;
            uweights(11,6 ,i)=   0.000001439560393714;
            uweights(11,7 ,i)=   0.4293597523561;
            uweights(11,8 ,i)=   0.1172278751677;
            uweights(11,9 ,i)=   0.01191139544491;
            uweights(11,10,i)=   0.0003468194663233;
            uweights(11,11,i)=   0.000001439560393714;
%
            unodes  (12,1 ,i)=   0.314240376254359;
            unodes  (12,2 ,i)=   0.947788391240164;
            unodes  (12,3 ,i)=   1.597682635152605;
            unodes  (12,4 ,i)=   2.279507080501060;
            unodes  (12,5 ,i)=   3.020637025120890;
            unodes  (12,6 ,i)=   3.889724897869782;
            unodes  (12,7 ,i)= - 0.314240376254359;
            unodes  (12,8 ,i)= - 0.947788391240164;
            unodes  (12,9 ,i)= - 1.597682635152605;
            unodes  (12,10,i)= - 2.279507080501060;
            unodes  (12,11,i)= - 3.020637025120890;
            unodes  (12,12,i)= - 3.889724897869782;
            uweights(12,1 ,i)=   0.5701352362625;
            uweights(12,2 ,i)=   0.2604923102642;
            uweights(12,3 ,i)=   0.05160798561588;
            uweights(12,4 ,i)=   0.003905390584629;
            uweights(12,5 ,i)=   0.00008573687043588;
            uweights(12,6 ,i)=   0.0000002658551684356;
            uweights(12,7 ,i)=   0.5701352362625;
            uweights(12,8 ,i)=   0.2604923102642;
            uweights(12,9 ,i)=   0.05160798561588;
            uweights(12,10,i)=   0.003905390584629;
            uweights(12,11,i)=   0.00008573687043588;
            uweights(12,12,i)=   0.0000002658551684356;
%
            unodes  (13,1, i)=   0;
            unodes  (13,2, i)=   0.605763879171060;
            unodes  (13,3, i)=   1.220055036590748;
            unodes  (13,4, i)=   1.853107651601512;
            unodes  (13,5, i)=   2.519735685678238;
            unodes  (13,6, i)=   3.246608978372410;
            unodes  (13,7, i)=   4.101337596178640;
            unodes  (13,8, i)= - 0.605763879171060;
            unodes  (13,9, i)= - 1.220055036590748;
            unodes  (13,10,i)= - 1.853107651601512;
            unodes  (13,11,i)= - 2.519735685678238;
            unodes  (13,12,i)= - 3.246608978372410;
            unodes  (13,13,i)= - 4.101337596178640;
            uweights(13,1, i)=   0.6043931879211;
            uweights(13,2, i)=   0.4216162968985;
            uweights(13,3, i)=   0.1403233206870;
            uweights(13,4, i)=   0.02086277529617;
            uweights(13,5, i)=   0.001207459992719;
            uweights(13,6, i)=   0.00002043036040271;
            uweights(13,7, i)=   0.00000004825731850073;
            uweights(13,8, i)=   0.4216162968985;
            uweights(13,9, i)=   0.1403233206870;
            uweights(13,10,i)=   0.02086277529617;
            uweights(13,11,i)=   0.001207459992719;
            uweights(13,12,i)=   0.00002043036040271;
            uweights(13,13,i)=   0.00000004825731850073;
%
            unodes  (14,1, i)=   0.29174551067256;
            unodes  (14,2, i)=   0.87871378732940;
            unodes  (14,3, i)=   1.47668273114114;
            unodes  (14,4, i)=   2.09518325850772;
            unodes  (14,5, i)=   2.74847072498540;
            unodes  (14,6, i)=   3.46265693360227;
            unodes  (14,7, i)=   4.30444857047363;
            unodes  (14,8, i)= - 0.29174551067256;
            unodes  (14,9, i)= - 0.87871378732940;
            unodes  (14,10,i)= - 1.47668273114114;
            unodes  (14,11,i)= - 2.09518325850772;
            unodes  (14,12,i)= - 2.74847072498540;
            unodes  (14,13,i)= - 3.46265693360227;
            unodes  (14,14,i)= - 4.30444857047363;
            uweights(14,1 ,i)=   0.5364059097121;
            uweights(14,2 ,i)=   0.2731056090642;
            uweights(14,3 ,i)=   0.06850553422347;
            uweights(14,4 ,i)=   0.007850054726458;
            uweights(14,5 ,i)=   0.0003550926135519;
            uweights(14,6 ,i)=   0.000004716484355019;
            uweights(14,7 ,i)=   0.000000008628591168125;
            uweights(14,8 ,i)=   0.5364059097121;
            uweights(14,9 ,i)=   0.2731056090642;
            uweights(14,10,i)=   0.06850553422347;
            uweights(14,11,i)=   0.007850054726458;
            uweights(14,12,i)=   0.0003550926135519;
            uweights(14,13,i)=   0.000004716484355019;
            uweights(14,14,i)=   0.000000008628591168125;
%
            unodes  (15,1, i)=   0;
            unodes  (15,2, i)=   0.56506958325558;
            unodes  (15,3, i)=   1.13611558521092;
            unodes  (15,4, i)=   1.71999257518649;
            unodes  (15,5, i)=   2.32573248617386;
            unodes  (15,6, i)=   2.96716692790560;
            unodes  (15,7, i)=   3.66995037340445;
            unodes  (15,8, i)=   4.49999070730939;
            unodes  (15,9, i)= - 0.56506958325558;
            unodes  (15,10,i)= - 1.13611558521092;
            unodes  (15,11,i)= - 1.71999257518649;
            unodes  (15,12,i)= - 2.32573248617386;
            unodes  (15,13,i)= - 2.96716692790560;
            unodes  (15,14,i)= - 3.66995037340445;
            unodes  (15,15,i)= - 4.49999070730939;
            uweights(15,1 ,i)=   0.5641003087264;
            uweights(15,2 ,i)=   0.4120286874989;
            uweights(15,3 ,i)=   0.1584889157959;
            uweights(15,4 ,i)=   0.03078003387255;
            uweights(15,5 ,i)=   0.002778068842913;
            uweights(15,6 ,i)=   0.0001000044412325;
            uweights(15,7 ,i)=   0.000001059115547711;
            uweights(15,8 ,i)=   0.00000000152247504254;
            uweights(15,9 ,i)=   0.4120286874989;
            uweights(15,10,i)=   0.1584889157959;
            uweights(15,11,i)=   0.03078003387255;
            uweights(15,12,i)=   0.002778068842913;
            uweights(15,13,i)=   0.0001000044412325;
            uweights(15,14,i)=   0.000001059115547711;
            uweights(15,15,i)=   0.00000000152247504254;
%
            unodes  (16,1, i)=   0.27348104613815;
            unodes  (16,2, i)=   0.82295144914466;
            unodes  (16,3, i)=   1.38025853919888;
            unodes  (16,4, i)=   1.95178799091625;
            unodes  (16,5, i)=   2.54620215784748;
            unodes  (16,6, i)=   3.17699916197996;
            unodes  (16,7, i)=   3.86944790486012;
            unodes  (16,8, i)=   4.68873893930582;
            unodes  (16,9, i)= - 0.27348104613815;
            unodes  (16,10,i)= - 0.82295144914466;
            unodes  (16,11,i)= - 1.38025853919888;
            unodes  (16,12,i)= - 1.95178799091625;
            unodes  (16,13,i)= - 2.54620215784748;
            unodes  (16,14,i)= - 3.17699916197996;
            unodes  (16,15,i)= - 3.86944790486012;
            unodes  (16,16,i)= - 4.68873893930582;
            uweights(16,1, i)=   0.5079294790166;
            uweights(16,2, i)=   0.2806474585285;
            uweights(16,3, i)=   0.08381004139899;
            uweights(16,4, i)=   0.01288031153551;
            uweights(16,5, i)=   0.0009322840086242;
            uweights(16,6, i)=   0.00002711860092538;
            uweights(16,7, i)=   0.0000002320980844865;
            uweights(16,8, i)=   0.0000000002654807474011;
            uweights(16,9, i)=   0.5079294790166;
            uweights(16,10,i)=   0.2806474585285;
            uweights(16,11,i)=   0.08381004139899;
            uweights(16,12,i)=   0.01288031153551;
            uweights(16,13,i)=   0.0009322840086242;
            uweights(16,14,i)=   0.00002711860092538;
            uweights(16,15,i)=   0.0000002320980844865;
            uweights(16,16,i)=   0.0000000002654807474011;
%
            unodes  (17,1 ,i)=   0;
            unodes  (17,2 ,i)=   0.5316330013427;
            unodes  (17,3 ,i)=   1.0676487257435;
            unodes  (17,4 ,i)=   1.6129243142212;
            unodes  (17,5 ,i)=   2.1735028266666;
            unodes  (17,6 ,i)=   2.7577629157039;
            unodes  (17,7 ,i)=   3.3789320911415;
            unodes  (17,8 ,i)=   4.0619466758755;
            unodes  (17,9 ,i)=   4.8713451936744;
            unodes  (17,10,i)= - 0.5316330013427;
            unodes  (17,11,i)= - 1.0676487257435;
            unodes  (17,12,i)= - 1.6129243142212;
            unodes  (17,13,i)= - 2.1735028266666;
            unodes  (17,14,i)= - 2.7577629157039;
            unodes  (17,15,i)= - 3.3789320911415;
            unodes  (17,16,i)= - 4.0619466758755;
            unodes  (17,17,i)= - 4.8713451936744;
            uweights(17,1 ,i)=   0.5309179376249;
            uweights(17,2 ,i)=   0.4018264694704;
            uweights(17,3 ,i)=   0.1726482976701;
            uweights(17,4 ,i)=   0.04092003414976;
            uweights(17,5 ,i)=   0.005067349957628;
            uweights(17,6 ,i)=   0.0002986432866978;
            uweights(17,7 ,i)=   0.000007112289140021;
            uweights(17,8 ,i)=   0.00000004977078981631;
            uweights(17,9 ,i)=   0.00000000004580578930799;
            uweights(17,10,i)=   0.4018264694704;
            uweights(17,11,i)=   0.1726482976701;
            uweights(17,12,i)=   0.04092003414976;
            uweights(17,13,i)=   0.005067349957628;
            uweights(17,14,i)=   0.0002986432866978;
            uweights(17,15,i)=   0.000007112289140021;
            uweights(17,16,i)=   0.00000004977078981631;
            uweights(17,17,i)=   0.00000000004580578930799;
%
            unodes  (18,1 ,i)=   0.2582677505191;
            unodes  (18,2 ,i)=   0.7766829192674;
            unodes  (18,3 ,i)=   1.3009208583896;
            unodes  (18,4 ,i)=   1.8355316042616;
            unodes  (18,5 ,i)=   2.3862990891667;
            unodes  (18,6 ,i)=   2.9613775055316;
            unodes  (18,7 ,i)=   3.5737690684863;
            unodes  (18,8 ,i)=   4.2481178735681;
            unodes  (18,9 ,i)=   5.0483640088745;
            unodes  (18,10,i)= - 0.2582677505191;
            unodes  (18,11,i)= - 0.7766829192674;
            unodes  (18,12,i)= - 1.3009208583896;
            unodes  (18,13,i)= - 1.8355316042616;
            unodes  (18,14,i)= - 2.3862990891667;
            unodes  (18,15,i)= - 2.9613775055316;
            unodes  (18,16,i)= - 3.5737690684863;
            unodes  (18,17,i)= - 4.2481178735681;
            unodes  (18,18,i)= - 5.0483640088745;
            uweights(18,1 ,i)=   0.4834956947255;
            uweights(18,2 ,i)=   0.2848072856700;
            uweights(18,3 ,i)=   0.09730174764132;
            uweights(18,4 ,i)=   0.01864004238754;
            uweights(18,5 ,i)=   0.001888522630268;
            uweights(18,6 ,i)=   0.00009181126867929;
            uweights(18,7 ,i)=   0.000001810654481093;
            uweights(18,8 ,i)=   0.00000001046720579579;
            uweights(18,9 ,i)=   0.000000000007828199772116;
            uweights(18,10,i)=   0.4834956947255;
            uweights(18,11,i)=   0.2848072856700;
            uweights(18,12,i)=   0.09730174764132;
            uweights(18,13,i)=   0.01864004238754;
            uweights(18,14,i)=   0.001888522630268;
            uweights(18,15,i)=   0.00009181126867929;
            uweights(18,16,i)=   0.000001810654481093;
            uweights(18,17,i)=   0.00000001046720579579;
            uweights(18,18,i)=   0.000000000007828199772116;
%
            unodes  (19,1 ,i)=   0;
            unodes  (19,2 ,i)=   0.5035201634239;
            unodes  (19,3 ,i)=   1.0103683871343;
            unodes  (19,4 ,i)=   1.5241706193935;
            unodes  (19,5 ,i)=   2.0492317098506;
            unodes  (19,6 ,i)=   2.5911337897945;
            unodes  (19,7 ,i)=   3.1578488183476;
            unodes  (19,8 ,i)=   3.7621873519640;
            unodes  (19,9 ,i)=   4.4285328066038;
            unodes  (19,10,i)=   5.2202716905375;
            unodes  (19,11,i)= - 0.5035201634239;
            unodes  (19,12,i)= - 1.0103683871343;
            unodes  (19,13,i)= - 1.5241706193935;
            unodes  (19,14,i)= - 2.0492317098506;
            unodes  (19,15,i)= - 2.5911337897945;
            unodes  (19,16,i)= - 3.1578488183476;
            unodes  (19,17,i)= - 3.7621873519640;
            unodes  (19,18,i)= - 4.4285328066038;
            unodes  (19,19,i)= - 5.2202716905375;
            uweights(19,1 ,i)=   0.5029748882762;
            uweights(19,2 ,i)=   0.3916089886130;
            uweights(19,3 ,i)=   0.1836327013070;
            uweights(19,4 ,i)=   0.05081038690905;
            uweights(19,5 ,i)=   0.007988866777723;
            uweights(19,6 ,i)=   0.0006708775214072;
            uweights(19,7 ,i)=   0.00002720919776316;
            uweights(19,8 ,i)=   0.0000004488243147223;
            uweights(19,9 ,i)=   0.000000002163051009864;
            uweights(19,10,i)=   0.000000000001326297094499;
            uweights(19,11,i)=   0.3916089886130;
            uweights(19,12,i)=   0.1836327013070;
            uweights(19,13,i)=   0.05081038690905;
            uweights(19,14,i)=   0.007988866777723;
            uweights(19,15,i)=   0.0006708775214072;
            uweights(19,16,i)=   0.00002720919776316;
            uweights(19,17,i)=   0.0000004488243147223;
            uweights(19,18,i)=   0.000000002163051009864;
            uweights(19,19,i)=   0.000000000001326297094499;
%
            unodes  (20,1 ,i)=   0.2453407083009;
            unodes  (20,2 ,i)=   0.7374737285454;
            unodes  (20,3 ,i)=   1.2340762153953; 
            unodes  (20,4 ,i)=   1.7385377121166;
            unodes  (20,5 ,i)=   2.2549740020893;
            unodes  (20,6 ,i)=   2.7888060584281;
            unodes  (20,7 ,i)=   3.3478545673832;
            unodes  (20,8 ,i)=   3.9447640401156;
            unodes  (20,9 ,i)=   4.6036824495507;
            unodes  (20,10,i)=   5.3874808900112;
            unodes  (20,11,i)= - 0.2453407083009;
            unodes  (20,12,i)= - 0.7374737285454;
            unodes  (20,13,i)= - 1.2340762153953;
            unodes  (20,14,i)= - 1.7385377121166;
            unodes  (20,15,i)= - 2.2549740020893;
            unodes  (20,16,i)= - 2.7888060584281;
            unodes  (20,17,i)= - 3.3478545673832;
            unodes  (20,18,i)= - 3.9447640401156;
            unodes  (20,19,i)= - 4.6036824495507;
            unodes  (20,20,i)= - 5.3874808900112;
            uweights(20,1 ,i)=   0.4622436696006;
            uweights(20,2 ,i)=   0.2866755053628;
            uweights(20,3 ,i)=   0.1090172060200;
            uweights(20,4 ,i)=   0.02481052088746;
            uweights(20,5 ,i)=   0.003243773342238;
            uweights(20,6 ,i)=   0.0002283386360163;
            uweights(20,7 ,i)=   0.000007802556478532;
            uweights(20,8 ,i)=   0.0000001086069370769;
            uweights(20,9 ,i)=   0.0000000004399340992273;
            uweights(20,10,i)=   0.0000000000002229393645534;
            uweights(20,11,i)=   0.4622436696006;
            uweights(20,12,i)=   0.2866755053628;
            uweights(20,13,i)=   0.1090172060200;
            uweights(20,14,i)=   0.02481052088746;
            uweights(20,15,i)=   0.003243773342238;
            uweights(20,16,i)=   0.0002283386360163;
            uweights(20,17,i)=   0.000007802556478532;
            uweights(20,18,i)=   0.0000001086069370769;
            uweights(20,19,i)=   0.0000000004399340992273;
            uweights(20,20,i)=   0.0000000000002229393645534;
%
            s2    =sqrt(2);
            spi   =sqrt(2*4*atan(1));
            s2_spi=s2/spi;
            for k1=1:mxhq
            for k2=1:k1
                unodes(k1,k2,i)=unodes(k1,k2,i)*s2;
                uweights(k1,k2,i)=uweights(k1,k2,i)*s2_spi;
            end
            end
%
         elseif iPDF(i)== 2 %Legendre
%
            unodes  (1,1,i)=   0;
            uweights(1,1,i)=   1;
%
            unodes  (2,1,i)=   1/sqrt(3);
            unodes  (2,2,i)= - 1/sqrt(3);
            uweights(2,1,i)=   0.5;
            uweights(2,2,i)=   0.5;
%
            unodes  (3,1,i)=   0;
            unodes  (3,2,i)=   sqrt(0.6);
            unodes  (3,3,i)= - sqrt(0.6);
            uweights(3,1,i)=   8/9 / 2;
            uweights(3,2,i)=   5/9 / 2;
            uweights(3,3,i)=   5/9 / 2;
%
            unodes  (4,1,i)= - 0.3399810435848563;
            unodes  (4,2,i)=   0.3399810435848563;
            unodes  (4,3,i)= - 0.8611363115940526;
            unodes  (4,4,i)=   0.8611363115940526;
            uweights(4,1,i)=   0.6521451548625461 /2;
            uweights(4,2,i)=   0.6521451548625461 /2;
            uweights(4,3,i)=   0.3478548451374538 /2;
            uweights(4,4,i)=   0.3478548451374538 /2;
%
            unodes  (5,1,i)=   0.0000000000000000;
            unodes  (5,2,i)= - 0.5384693101056831;
            unodes  (5,3,i)=   0.5384693101056831;
            unodes  (5,4,i)= - 0.9061798459386640;
            unodes  (5,5,i)=   0.9061798459386640;
            uweights(5,1,i)=   0.5688888888888889 /2;
            uweights(5,2,i)=   0.4786286704993665 /2;
            uweights(5,3,i)=   0.4786286704993665 /2;
            uweights(5,4,i)=   0.2369268850561891 /2;
            uweights(5,5,i)=   0.2369268850561891 /2;
%
            unodes  (6,1,i)=   0.6612093864662645;
            unodes  (6,2,i)= - 0.6612093864662645;
            unodes  (6,3,i)= - 0.2386191860831969;
            unodes  (6,4,i)=   0.2386191860831969;
            unodes  (6,5,i)= - 0.9324695142031521;
            unodes  (6,6,i)=   0.9324695142031521;
            uweights(6,1,i)=   0.3607615730481386 /2;
            uweights(6,2,i)=   0.3607615730481386 /2;
            uweights(6,3,i)=   0.4679139345726910 /2;
            uweights(6,4,i)=   0.4679139345726910 /2;
            uweights(6,5,i)=   0.1713244923791704 /2;
            uweights(6,6,i)=   0.1713244923791704 /2;
%
            unodes  (7,1,i)=   0;
            unodes  (7,2,i)=   0.4058451513773972;
            unodes  (7,3,i)= - 0.4058451513773972;
            unodes  (7,4,i)= - 0.7415311855993945;
            unodes  (7,5,i)=   0.7415311855993945;
            unodes  (7,6,i)=   0.9491079123427585;
            unodes  (7,7,i)= - 0.9491079123427585;
            uweights(7,1,i)=   0.4179591836734694 /2;
            uweights(7,2,i)=   0.3818300505051189 /2;
            uweights(7,3,i)=   0.3818300505051189 /2;
            uweights(7,4,i)=   0.2797053914892766 /2;
            uweights(7,5,i)=   0.2797053914892766 /2;
            uweights(7,6,i)=   0.1294849661688697 /2;
            uweights(7,7,i)=   0.1294849661688697 /2;
%
            unodes(8,1,i)  =  -0.1834346424956498;
            unodes(8,2,i)  =   0.1834346424956498;
            unodes(8,3,i)  =  -0.5255324099163290;
            unodes(8,4,i)  =   0.5255324099163290;
            unodes(8,5,i)  =  -0.7966664774136267;
            unodes(8,6,i)  =   0.7966664774136267;
            unodes(8,7,i)  =  -0.9602898564975363;
            unodes(8,8,i)  =   0.9602898564975363;
            uweights(8,1,i)=   0.3626837833783620/2;
            uweights(8,2,i)=   0.3626837833783620/2;
            uweights(8,3,i)=   0.3137066458778873/2;
            uweights(8,4,i)=   0.3137066458778873/2;
            uweights(8,5,i)=   0.2223810344533745/2;
            uweights(8,6,i)=   0.2223810344533745/2;
            uweights(8,7,i)=   0.1012285362903763/2;
            uweights(8,8,i)=   0.1012285362903763/2;
%
            unodes(9,1,i)  =   0;
            unodes(9,2,i)  =  -0.8360311073266358;
            unodes(9,3,i)  =   0.8360311073266358;
            unodes(9,4,i)  =  -0.9681602395076261;
            unodes(9,5,i)  =   0.9681602395076261;
            unodes(9,6,i)  =  -0.3242534234038089;
            unodes(9,7,i)  =   0.3242534234038089;
            unodes(9,8,i)  =  -0.6133714327005904;
            unodes(9,9,i)  =   0.6133714327005904;
            uweights(9,1,i)=   0.3302393550012598/2;
            uweights(9,2,i)=   0.1806481606948574/2;
            uweights(9,3,i)=   0.1806481606948574/2;
            uweights(9,4,i)=   0.0812743883615744/2;
            uweights(9,5,i)=   0.0812743883615744/2;
            uweights(9,6,i)=   0.3123470770400029/2;
            uweights(9,7,i)=   0.3123470770400029/2;
            uweights(9,8,i)=   0.2606106964029354/2;
            uweights(9,9,i)=   0.2606106964029354/2;
%
            unodes(10,1,i)   =  -0.1488743389816312;
            unodes(10,2,i)   =   0.1488743389816312;
            unodes(10,3,i)   =  -0.4333953941292472;
            unodes(10,4,i)   =   0.4333953941292472;
            unodes(10,5,i)   =  -0.6794095682990244;
            unodes(10,6,i)   =   0.6794095682990244;
            unodes(10,7,i)   =  -0.8650633666889845;
            unodes(10,8,i)   =   0.8650633666889845;
            unodes(10,9,i)   =  -0.9739065285171717;
            unodes(10,10,i)  =   0.9739065285171717;
            uweights(10,1,i) =   0.2955242247147529/2;
            uweights(10,2,i) =   0.2955242247147529/2;
            uweights(10,3,i) =   0.2692667193099963/2;
            uweights(10,4,i) =   0.2692667193099963/2;
            uweights(10,5,i) =   0.2190863625159820/2;
            uweights(10,6,i) =   0.2190863625159820/2;
            uweights(10,7,i) =   0.1494513491505806/2;
            uweights(10,8,i) =   0.1494513491505806/2;
            uweights(10,9,i) =   0.0666713443086881/2;
            uweights(10,10,i)=   0.0666713443086881/2;
%
         elseif iPDF(i)== 3  %Laguerre
%
            unodes  (1,1,i)=   1;
            uweights(1,1,i)=   1;
%
            unodes  (2,1,i)=   0.5857864376269051;
            unodes  (2,2,i)=   3.414213562373095;
            uweights(2,1,i)=   0.8535533905932737;
            uweights(2,2,i)=   0.1464466094067262;
%
            unodes  (3,1,i)=   0.4157745567834789;
            unodes  (3,2,i)=   2.294280360279041;
            unodes  (3,3,i)=   6.289945082937480;
            uweights(3,1,i)=   0.7110930099291729;
            uweights(3,2,i)=   0.2785177335692409;
            uweights(3,3,i)=   0.1038925650158615e-1;
%
            unodes  (4,1,i)=   0.3225476896193923;
            unodes  (4,2,i)=   1.745761101158346;
            unodes  (4,3,i)=   4.536620296921128;
            unodes  (4,4,i)=   9.395070912301136;
            uweights(4,1,i)=   0.6031541043416337;
            uweights(4,2,i)=   0.3574186924377999;
            uweights(4,3,i)=   0.3888790851500538e-1;
            uweights(4,4,i)=   0.5392947055613278e-3;
%
            unodes  (5,1,i)=   0.2635603197181410;
            unodes  (5,2,i)=   1.413403059106517;
            unodes  (5,3,i)=   3.596425771040721;
            unodes  (5,4,i)=   7.085810005858835;
            unodes  (5,5,i)=   12.64080084427578;
            uweights(5,1,i)=   0.5217556105828089;
            uweights(5,2,i)=   0.3986668110831760;
            uweights(5,3,i)=   0.7594244968170749e-1;
            uweights(5,4,i)=   0.3611758679922050e-2;
            uweights(5,5,i)=   0.2336997238577621e-4;
%
        end
%     =====
      end
%     =====
%???????????????????????????????????????????????????????????????????????
if mxs<((ko+1)^id)
    disp('Full Tensorized Grid Used. Should be Smolyak.')
    n    = ko+1;
    mxs  = (ko+1)^id;
    snodes = zeros(mxs,id);
    sweights = zeros(1,mxs); 
    % do I need preallocation for psi = zeros(nhe+1,mxs);
    for ihq = 1:mxs
        sweights(ihq) = 1;
        for i = 1:id
            if i == 1
                itensor = mod(ihq-1,n)+1;
            else
                itensor = mod(floor((ihq-1)/(n^(i-1))),n)+1;
            end 
            snodes(ihq,i) = unodes(n,itensor,i);
            sweights(ihq) = sweights(ihq)*uweights(n,itensor,i);
        end
    end %for
else
    disp('Full Tensorized Grid Used.')
    n    = ko+1;
    mxs  = (ko+1)^id;
    snodes = zeros(mxs,id);
    sweights = zeros(1,mxs); 
    % do I need preallocation for psi = zeros(nhe+1,mxs);
    for ihq = 1:mxs
        sweights(ihq) = 1;
        for i = 1:id
            if i == 1
                itensor = mod(ihq-1,n)+1;
            else
                itensor = mod(floor((ihq-1)/(n^(i-1))),n)+1;
            end 
            snodes(ihq,i) = unodes(n,itensor,i);
            sweights(ihq) = sweights(ihq)*uweights(n,itensor,i);
        end
    end %for
end %mxs<(ko+1)^d
%????????????????????????????????????????????????????????????????????????
[psi,nbox] = Hermite(ko,snodes,mxs,id,nhe,iPDF);
end

function [psi,nbox] = Hermite(ko,snodes,mxs,id,nhe,iPDF)
%calculates the values of all id-dimensional polynomials
%at all quadrature nodes (snodes) and stores them at psi
%???????????????????????????????????????????????????????????????????????
nbox=zeros(nhe+1,id);
kn=0;
for ik=0:ko
    mxind=nchoosek(ik+id-1,ik);
    %ind=zeros(mxind,id);
    [kount,ind] = drop(id,ik+id,mxind);
    for i=1:kount
        kn=kn+1;
        for j=1:id
            nbox(kn,j)=ind(i,j)-1;
        end    
    end
end

%???????????????????????????????????????????????????????????????????????
psi = zeros(nhe+1,mxs);
for i1=0:nhe
for i2=1:mxs 
    prod=1;
    for i3=1:id
        prod=prod*He(nbox(i1+1,i3),snodes(i2,i3),iPDF(i3));
    end
    psi(i1+1,i2)=prod;
end
end
end %function

%=======================
function He = He(n,x,kPDF) %Unitary orthogonal polynomials
%=======================
if kPDF == 1 %propabilistic Hermite
   if n==0
       He=1;
       return
   elseif n==1
       He=x;
       return
   end
   h0=1; h1=x;
   for i=2:n
       h=x*h1-(i-1)*h0;
       h0=h1;
       h1=h;
   end
   He = h/sqrt(factorial(n));   
elseif kPDF==2 %Legendre
    if n==0
        He=1;
        return
    elseif n==1
        He=x*sqrt(3);
       return 
    end
    h0=1; h1=x; %do we need an sqrt(3) here?
    for i=2:n
        j=i-1;
        h=(2*j+1)*x*h1-j*h0;
        h0=h1;
        h1=h;
    end
    anorm=1/(2*n+1);
    He=h/sqrt(anorm);    
elseif kPDF==3 %Laguerre
    if n==0 
        He=1;
        return
    elseif n==1
        He=1-x;
        return
    end
    h0=1; h1=1-x;
    for i=2:n
        ic=i-1;
        c1=(2*ic+1)-x;      %Where do we use that?
        c0=c1*h1+c0*h0;     %Where do we use that?
        h=h/dfloat(ic+1);   %Does not look like (2.10)
        h0=h1;
        h1=h;
    end
    He=h;
end   
end

%???????????????????????????????????????????????????????????????????????
function [kount,ind] = drop(id,l,mxind)
%Thomas Derstner drop algorithm
%Generator of multi-indices a=(a(1),a(2),...,a(id))
%such that a(i)>=1 for i=1,2,...,d and a(1)+a(2)+...a(id)=l
%inputs : id,l
%output : kount (# of multi-indices), ind(kount,id)
ind=zeros(mxind,id);
k  =zeros(id); 
kh =zeros(id);
ip=1; kount=0;
for i=1:id
    k(i)=1;
    kh(i)=l-id+1;
end
if l<id 
    return
end
if l == id
    kount=1;
    for i=1:id
        ind(kount,i)=1;
    end
    return
end
while (k(id)<=l)
    k(ip)=k(ip)+1; 
    if (k(ip)>kh(ip))
       if (ip~=id)
           k(ip)=1; ip=ip+1;
       end
    else
       for i=1:ip-1
           kh(i)=kh(ip)-k(ip)+1;
       end %for
       k(1)=kh(1);
       ip=1;
       kount=kount+1;
       for i=1:id
           ind(kount,i)=k(i);
       end %for
    end %if        
end    %while
ii=nchoosek(l-1,id-1);
if kount~=ii 
    disp('Error in subroutine drop')
    return
end %if
return 
end
%???????????????????????????????????????????????????????????????????????    

