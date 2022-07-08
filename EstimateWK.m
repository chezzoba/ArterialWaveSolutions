addpath("Algorithms/");

RW1s = [0,0,0,0,0,0,0,0,0,0,5.1918,19.1515,9.882,11.7617,17.4352,34.1378,34.1378,74.0167,5.9149,5.9149]*(10^7);
RW2s = [0,0,0,0,0,0,0,0,0,0,10.6080,52.2129,13.0183,7.5726,5.5097,5.3949,5.3949,46.2252,10.1737,10.1737]*(10^8);
CWKs = [0,0,0,0,0,0,0,0,0,0,8.6974,1.767,7.0871,12.1836,16.7453,17.1017,17.1017,1.9959,9.0686,9.0686]*(10^-10);

number = 11;
flows = load(sprintf('PreviousCode/automeris/aorta/bc_outlet_%d_flow.csv', number));
pressures = load(sprintf('PreviousCode/automeris/aorta/bc_outlet_%d_Pressure.csv', number));

xpred = WindKessels(pressures, flows);

xact = [RW1s(11), RW2s(11), CWKs(11)]

Perr = 100.*ConOptimisation.relerr(xpred, xact)