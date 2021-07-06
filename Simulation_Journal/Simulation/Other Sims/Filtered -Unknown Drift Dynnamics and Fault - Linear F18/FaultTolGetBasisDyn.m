function [phiF,g,fo] = FaultTolGetBasisDyn(x)

phiF = [x(1);x(2)]; %drift dynamic basis
Wg = [-0.194,-0.03593;
       -19.29,-3.803]'; %control effectiveness actual weight
phiG = eye(2); %control effectiveness basis

g = Wg'*phiG;

fo = zeros(size(x,1),1); %Known part of dynamics



end