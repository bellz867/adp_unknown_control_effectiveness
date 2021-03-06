function [pts] = FaultTolADPExtTraj(t,x,auxdata)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
m = length(auxdata.Rvec);
n = length(x);
nmx = ((x'*x)+auxdata.dbar2)/(1+auxdata.scale2*(x'*x)); %Skrinking factor
% nmx=1;
% pts = unifrnd(-nmx,nmx,n,auxdata.numpoints);
% pts = auxdata.BEscale*nmx*randn(n,auxdata.numpoints);


pts = zeros(size(x,1),auxdata.numpoints);

for i=1:length(x)
    for j=1:auxdata.numpoints
    pts(i,j) = auxdata.BEscale*nmx*sum(sin(j*auxdata.freq(i,:)*t+auxdata.phase(i,:)));
%     pts(i,j) = unifrnd(-auxdata.phase(i,j)*auxdata.BEscale*nmx,auxdata.phase(i,j)*auxdata.BEscale*nmx);
    end
end
    
    


%     XC1=sum(sin(freq(1,:)*t+Phase(1,:)));
%     XC2=sum(sin(freq(2,:)*t+Phase(2,:)));
%     XC3=sum(sin(freq(3,:)*t+Phase(3,:)));
%     XC4=sum(sin(freq(4,:)*t+Phase(4,:)));
%     XC1=((1/NN)*sin(6*pi*t)...
%         +(1/NN)*sin(exp(5)*t)...
%         +(1/NN)*cos(10*sqrt(pi)*t)...
%         +(1/NN)*sin(sqrt(2)*t)...
%         +(1/NN)*cos(5*exp(pi)*t)...
%         +(1/NN)*sin(pi^(exp(1))*t)...
%         +(1/NN)*cos(2*t)...
%         +(1/NN)*sin(sqrt(3)*t)...
%         +(1/NN)*cos(2*sqrt(exp(1))*t));
%     XC2=((1/NN)*cos(2.5*pi*t/2)...
%         +(1/NN)*sin(exp(0.5)*t)...
%         +(1/NN)*cos(2.5*sqrt(1.5*pi)*t)...
%         +(1/NN)*sin(sqrt(5/2)*t)...
%         +(1/NN)*cos(2.5*exp(sqrt(pi))*t)...
%         +(1/NN)*sin(2*pi^(exp(2))*t)...
%         +(1/NN)*sin(t)...
%         +(1/NN)*cos(sqrt(3.5)*t)...
%         +(1/NN)*sin(sqrt(exp(1.5))*t));
%     XC3=((1/NN)*sin(0.4*pi*t)...
%         +(1/NN)*cos(exp(3)*t)...
%         +(1/NN)*sin(2.5*sqrt(2*pi)*t)...
%         +(1/NN)*sin(sqrt(0.5)*t)...
%         +(1/NN)*cos(2.5*exp(pi/2)*t)...
%         +(1/NN)*sin(pi^(exp(1))*t)...
%         +(1/NN)*sin(2*t)...
%         +(1/NN)*cos(sqrt(7)*t)...
%         +(1/NN)*sin(2.5*((exp(1))^(1/3))*t));
%     XC4=((1/NN)*sin(pi*t)...
%         +(1/NN)*cos(2.5*exp(1)*t)...
%         +(1/NN)*sin((pi^(1/5))*t)...
%         +(1/NN)*sin(3*(2^(1/5))*t)...
%         +(1/NN)*sin(2.5*exp(pi^(1/3))*t)...
%         +(1/NN)*cos(pi^(exp(1*pi))*t)...
%         +(1/NN)*sin(2.5*t)...
%         +(1/NN)*sin(2.5*sqrt(3)*t)...
%         +(1/NN)*cos(sqrt(exp(1.5))*t));
%     excite=[XC1;XC2;XC3;XC4];
end

