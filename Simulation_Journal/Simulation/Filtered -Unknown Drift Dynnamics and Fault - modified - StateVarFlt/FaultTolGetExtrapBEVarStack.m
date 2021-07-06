function [criticSumNew,actorSumNew,ADPGammaSumNew,ADPStackSumNew] = FaultTolGetExtrapBEVarStack(time,x,thetaHat,WcHat,WaHat,GammaADP,muHat,auxdata)
%-------------------------------------------------------------------------%
%------------------ Build Bellman Extrapolation Library ------------------%
%-------------------------------------------------------------------------%

m = length(auxdata.Rvec);
n = length(x);
criticSum = zeros(length(WcHat),1);
actorSum = zeros(length(WcHat),1);
ADPGammaSum = zeros(length(WcHat),length(WcHat));
ADPStackSum = zeros(length(WcHat),length(WcHat));


nmx = ((x'*x)+auxdata.dbar2)/(1+auxdata.scale2*(x'*x)); %Skrinking factor
% nmx=1;
if auxdata.MovePoints == 0
    temp=linspace(-auxdata.BEscale*nmx,auxdata.BEscale*nmx,auxdata.numpoints);
    [temp1, temp2] = ndgrid(temp,temp);
    pts=[temp1(:) temp2(:)]'; % matrix with columns as selected points
    M = (auxdata.numpoints)^2; %Number of points
else
    pts = unifrnd(-auxdata.BEscale*nmx,auxdata.BEscale*nmx,n,auxdata.numpoints);
    % pts = auxdata.BEscale*nmx*randn(n,auxdata.numpoints);
    % pts = FaultTolADPExtTraj(time,x,auxdata);
    M =auxdata.numpoints;
end
% WfHat= thetaHat(1:n*auxdata.pf,:,:); %initialize drift dynamics weight estimate
%     WmuHat(:,:,t+1) = thetaHat(auxdata.pf+1:auxdata.pf+n*m,:,t+1); %initialize control effective + fault weight estimate
WgHat = reshape(thetaHat(n*auxdata.pf+1:auxdata.pGamma,:,:),auxdata.pg,m); 
     
for k = 1:auxdata.numpoints;
    xk = pts(:,k)+x; %Extrpolated state 
    %Get Basis 
    [phifk,phigk,gk,fok] = FaultTolGetBasisDyn(xk,auxdata); %Get basis functions for estimator 
    gHatk = gk*diag(WgHat'*phigk);
    [~,sigPrimek] = FaultTolGetBasisADP(xk,auxdata); %Get StaF basis functions and their gradients
    [ uk, uCostk,Gak,~,~] = FaultTolGetInput(xk,gHatk,WaHat,sigPrimek,auxdata);%Get input input cost and musat(sigPrime(sgn()-Tanh())) term for WaHat
    PHIk = [kron(phifk',eye(2)),kron(phigk',gk*diag(uk))]; %Get Combined CL basis
    Gk = 0; %Get combined g*u
    [~,omegak,~,deltak,rhok] = FaultTolGetADPVar(xk,fok,sigPrimek,uCostk,WcHat,thetaHat,PHIk,auxdata,2); %ADP estapolated parameters
    criticSum = criticSum + GammaADP*omegak*(deltak)/(rhok^2);
    actorSum = actorSum + Gak*omegak'*WcHat/(rhok^2);
    ADPGammaSum = ADPGammaSum + GammaADP*(omegak*omegak')/(rhok^2)*GammaADP;
    ADPStackSum = ADPStackSum + (omegak*omegak')/(rhok^2); 
end

criticSumNew = criticSum/auxdata.numpoints;
actorSumNew = actorSum/auxdata.numpoints;
ADPStackSumNew = ADPStackSum/auxdata.numpoints;
ADPGammaSumNew = ADPGammaSum/auxdata.numpoints;
 
end