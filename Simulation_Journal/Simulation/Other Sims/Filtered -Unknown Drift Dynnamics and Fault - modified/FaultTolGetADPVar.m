function   [fgHat,omega,instCost,delta,rho] = FaultTolGetADPVar(x,fo,sigPrime,Ucost,WcHat,thetaHat,PHI,auxdata,Type)

if Type==1
    nu = auxdata.nu1;
else
    nu = auxdata.nu2;
end

Sx = (1-auxdata.InputKp2*(x'*x))/(1+auxdata.InputKp2*(x'*x))^3*x'*auxdata.InputKp1;
fgHat = (PHI./auxdata.NormDyn)*thetaHat+(fo./auxdata.NormDyn); %approximate Dynamics
omega = sigPrime*fgHat; %regressor
instCost = x'*auxdata.Qx*x+Ucost; %instantaneous cost
delta = omega'*WcHat+instCost+Sx*fgHat; %Bellman Error
rho = (1+nu*(omega'*omega)); %normalization term
      
end
