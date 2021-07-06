function [phiF,g,fo] = FaultTolGetBasisDyn(x)

phiF = [x(1);x(2);x(2)*(1-(cos(2*x(1))+2)^2)]; %drift dynamic basis
Wg = [0,1;0,2]; %control effectiveness actual weight
phiG = [cos(2*x(1));1]; %control effectiveness basis

g = Wg'*phiG;

fo = zeros(size(x,1),1); %Known part of dynamics



end