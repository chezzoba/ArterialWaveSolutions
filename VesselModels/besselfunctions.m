%% Calculation of J and Y
function [J_13,Y_13,J_43,Y_43,f,z] = besselfunctions(a,s,omega,rho,beta)
%Firstly z needs to be calculated

%Defining the constant parameters of the problem

f = beta*(tan(a)^2)*sin(a)/(1-cos(a));    %eq. 2 (Papadakis 2019)
z = (2/3).*omega.*((rho.*f)^0.5)*s^(3/2);              %(Papadakis 2019)

J_13 = besselj(1/3,z);                             
Y_13 = bessely(1/3,z);   

J_43 = besselj(4/3,z);                             
Y_43 = bessely(4/3,z);                             

end

