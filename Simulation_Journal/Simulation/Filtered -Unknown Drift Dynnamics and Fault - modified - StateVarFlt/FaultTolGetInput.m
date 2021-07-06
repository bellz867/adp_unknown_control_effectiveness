function [ u, uCost,Ga,regInput,auxInput] = FaultTolGetInput(y,gin,WaHat,sigPrime,auxdata)
% this function calculates the input and cost of input
Rvec = auxdata.Rvec;
NormDyn = auxdata.NormDyn;
InputKp1 = auxdata.InputKp1; 
InputKp2 = auxdata.InputKp2;
m = length(Rvec);
sat = auxdata.sat; %input saturation
R = diag(Rvec);
% x = y/NormDyn; %Normalize states
x = y;
g = gin/NormDyn; %Normalize dynamics

if (sat == 0)
    
    addInput = (1-InputKp2*(y'*y))/(1+InputKp2*(y'*y))^3*x'*InputKp1;
    u = -0.5*(R\g')*(sigPrime'*WaHat+addInput');
    auxInput= -0.5*(R\g')*(addInput'); %Just the auziliary Input
    regInput = -0.5*(R\g')*(sigPrime'*WaHat);%Input without Auxiliary term

    uCost = u'*R*u;
    Ga = (sigPrime*g*(R\g')*sigPrime')'*WaHat/4;
else
    addInput = (1-InputKp2*(y'*y))/(1+InputKp2*(y'*y))^3*x'*InputKp1;
    us = -sat*tanh((1/sat)*(0.5*(R\g')*(sigPrime'*WaHat+addInput')));
    auxInput= -sat*tanh((1/sat)*0.5*(R\g')*(addInput')); %Just the auziliary Input
    regInput = -sat*tanh((1/sat)*0.5*(R\g')*(sigPrime'*WaHat));%Input without Auxiliary term

%     u = -sat*tanh((1/sat)*(0.5*(R\g')*(sigPrime'*WaHat)));
    u = us - ((sat - abs(us)) < eps(us)).*sign(us).*eps(us);
    
%     u(1) = u(1) - ...
%         ((sat - abs(u(1))) < eps(u(1))).*sign(u(1)).*eps(u(1));
%     u(2) = u(2) - ...
%         ((sat - abs(u(2))) < eps(u(2))).*sign(u(2)).*eps(u(2));
    
%     Ucost = 2*sat*R(1,1)*u(1)*atanh(u(1)/sat) + ...
%         sat^2*R(1,1)*log(1 - u(1)^2/sat^2)+...
%         2*sat*R(2,2)*u(2)*atanh(u(2)/sat) + ...
%         sat^2*R(2,2)*log(1 - u(2)^2/sat^2);
    uCost = 2*sat*u'*R*atanh(u/sat)+ sat^2*Rvec*log(ones(m)-u.^2./sat^2);
    Ga = sigPrime*(sat*g*tanh(1/(sat*0.02)*g'*(sigPrime'*WaHat+addInput'))+g*u);

end

