function [ u, uCost,Ga,Gc] = FaultTolGetInput(y,gin,WaHat,sigPrime,auxdata)
% this function calculates the input and cost of input
m = length(auxdata.Rvec);
sat = auxdata.sat; %input saturation
% x = y./auxdata.NormDyn; %Normalize states
x = y;
g = gin./auxdata.NormDyn; %Normalize dynamics

if (sat == 0)
    addInput = auxdata.InputKp1*(y'*y/(1+auxdata.InputKp2*(y'*y)))*x;
    u = -0.5*(auxdata.R\g')*(sigPrime'*WaHat+addInput);
%     u = -0.5*(auxdata.R\g')*(sigPrime'*WaHat);

    uCost = u'*auxdata.R*u;
    Ga = (sigPrime*g*(auxdata.R\g')*sigPrime')'*WaHat/4;
else
    addInput = auxdata.InputKp1*(y'*y/(1+auxdata.InputKp2*(y'*y)))*x;
    u = -sat*tanh((1/sat)*(0.5*(auxdata.R\g')*(sigPrime'*WaHat+addInput)));
%     u = -sat*tanh((1/sat)*(0.5*(auxdata.R\g')*(sigPrime'*WaHat)));
    u = u - ((sat - abs(u)) < eps(u)).*sign(u).*eps(u);
    
%     u(1) = u(1) - ...
%         ((sat - abs(u(1))) < eps(u(1))).*sign(u(1)).*eps(u(1));
%     u(2) = u(2) - ...
%         ((sat - abs(u(2))) < eps(u(2))).*sign(u(2)).*eps(u(2));
    
%     Ucost = 2*sat*auxdata.R(1,1)*u(1)*atanh(u(1)/sat) + ...
%         sat^2*auxdata.R(1,1)*log(1 - u(1)^2/sat^2)+...
%         2*sat*auxdata.R(2,2)*u(2)*atanh(u(2)/sat) + ...
%         sat^2*auxdata.R(2,2)*log(1 - u(2)^2/sat^2);
    uCost = 2*sat*u'*auxdata.R*atanh(u/sat)+ sat^2*auxdata.Rvec*log(ones(m)-u.^2./sat^2);

%     Ga =
%     -sigPrime*(sat*g*sign((1/(sat))*(0.5*(auxdata.R\g')*(sigPrime'*WaHat+auxdata.InputKp*x)))+g*u); %Usig signum
    Ga = -sigPrime*(sat*g*tanh(1/(sat*0.02)*g'*(sigPrime'*WaHat+addInput))+g*u);  %Using close approximation of signum
%     Ga = zeros(size(Ga));
    Gc = (u-sat*sign((1/sat)*(0.5*(auxdata.R\g')*(sigPrime'*WaHat+addInput))))'*g'*addInput;

end

