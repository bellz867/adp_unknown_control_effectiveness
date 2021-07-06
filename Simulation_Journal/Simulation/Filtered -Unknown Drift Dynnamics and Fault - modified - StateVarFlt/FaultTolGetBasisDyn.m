function [phiF,phiG,g,fo] = FaultTolGetBasisDyn(x,auxdata)

Wg = [0,1;0,2]; %control effectiveness actual weight
psiG = [cos(2*x(1));1]; %control effectiveness basis
g = Wg'*psiG;

if auxdata.ExactBasis==1
    phiF = [x(1);x(2);x(2)*(1-(cos(2*x(1))+2)^2)]; %drift dynamic basis
    
    if auxdata.ExactFaultBasis==1
        phiG = [1];
    else
%         phiG = [1,tanh(auxdata.Vg'*x)']';
%         phiG = [1,tanh(x)',cos(x)',sin(x)',exp(-x'*x)]';
%         phiG = [1,(log(1+exp(x)))']';
        phiG = [1,(log(1+exp(auxdata.Vg'*[1,x']')))']';
%         phiG = [1,cos(x)',sin(x)']';
        if auxdata.ComparePrev==1
            phiG = phiG(1,:);           
        end
    end

    
    fo = zeros(size(x,1),1); %Known part of dynamics
else
    
    xCenter = [-2,-2/3,0,2/3,2,-2,-2/3,0,2/3,2,-2,-2/3,0,2/3,2,-2,-2/3,0,2/3,2,-2,-2/3,0,2/3,2];
    yCenter = [-2,-2,-2,-2,-2,-2/3,-2/3,-2/3,-2/3,-2/3,0,0,0,0,0,2/3,2/3,2/3,2/3,2/3,2,2,2,2,2];
    xCenter = [0, 0,    0.8660,   -0.8660];
    yCenter = [0, 1.0000,   -0.5000,   -0.5000];
    
    Centers = [xCenter;yCenter];
    
    phiFGaus = zeros(1,length(yCenter));
    for i=1:length(yCenter)
        phiFGaus(i) = exp(-(x-Centers(:,i))'*(x-Centers(:,i)));
    end
    
    phiF = [1,tanh(x)',phiFGaus]';
    phiF = [1,phiFGaus]';
    
    
    if auxdata.ExactFaultBasis ==1
        phiG = [1];
    else
        phiG = [1,(log(1+exp(auxdata.Vg'*[1;x])))']';
    end
    
    fo = zeros(size(x,1),1); %Known part of dynamics
end



end