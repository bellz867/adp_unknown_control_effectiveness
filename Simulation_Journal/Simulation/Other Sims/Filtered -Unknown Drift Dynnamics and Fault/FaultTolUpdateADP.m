function [WcHatNew,WaHatNew,GammaNew] = FaultTolUpdateADP(WcHat,WaHat,Gamma,Ga,Gc,...
    criticSum,actorSum,ADPGammaSum,omega,delta,rho,auxdata,dt)



WcHatNew = WcHat + dt*DWcHat; %Update critic weight
WaHatNew = WaHat + dt*DWaHat; %Update actor weight
GammaNew = Gamma + dt*DGamma; %Update least squares gain

%Check thetaHat lies within bounds
[I1,J1] = find(WaHatNew < auxdata.WaHatLower);
[I2,J2] = find(WaHatNew > auxdata.WaHatUpper);

if isempty([I1,J1])~=1;
    for i=1:length(I1)
    WaHatNew (I1(i),J1(i)) = auxdata.WaHatLower;
    end
end

if isempty([I2,J2])~=1;
    for i=1:length(I2)
    WaHatNew (I2(i),J2(i)) = auxdata.WaHatUpper;
    end
end




end