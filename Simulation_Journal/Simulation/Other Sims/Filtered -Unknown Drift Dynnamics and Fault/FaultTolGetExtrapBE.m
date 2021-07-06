function [criticSum,actorSum,ADPStackSum] = FaultTolGetExtrapBE (time,x,thetaHat,WcHat,WaHat,muHat,auxdata)
%-------------------------------------------------------------------------%
%------------------ Build Bellman Extrapolation Library ------------------%
%-------------------------------------------------------------------------%

m = length(auxdata.Rvec);
n = length(x);
criticSum = zeros(length(WcHat),1);
actorSum = zeros(length(WcHat),length(WcHat));
ADPStackSum = zeros(length(WcHat),length(WcHat));

nmx = ((x'*x)+auxdata.dbar2)/(1+auxdata.scale2*(x'*x)); %Skrinking factor
% nmx=1;
% pts = unifrnd(-auxdata.BEscale*nmx,auxdata.BEscale*nmx,n,auxdata.numpoints);
pts = auxdata.BEscale*nmx*randn(n,auxdata.numpoints);

% pts = FaultTolADPExtTraj(time,x,auxdata);
    
for k = 1:auxdata.numpoints
    xk = pts(:,k)+x; %Extrpolated state  
    %Get Basis 
    [phifk,gok] = FaultTolGetBasisDyn(xk); %Get basis functions for estimator 
    gHatk = gok*diag(muHat);
    [~,sigPrimek] = FaultTolGetBasisADP(xk,auxdata); %Get StaF basis functions and their gradients
    [ uk, uCostk,Gak,Gck ] = FaultTolGetInput(gHatk,WaHat,sigPrimek,auxdata);%Get input input cost and musat(sigPrime(sgn()-Tanh())) term for WaHat
    PHIk = [kron(phifk',eye(2)),-gok*diag(uk)]; %Get Combined CL basis
    Gk = 0; %Get combined g*u
    [~,omegak,deltak,rhok] = FaultTolGetADPVar(xk,Gk,sigPrimek,uCostk,WcHat,thetaHat,PHIk,auxdata,2); %ADP estapolated parameters
    criticSum = criticSum + omegak*(deltak+Gck*WaHat)/rhok^2;
    actorSum = actorSum + Gak*omegak'/rhok^2;
    ADPStackSum = ADPStackSum + (omegak*omegak')/(rhok^(2));
end

 
end