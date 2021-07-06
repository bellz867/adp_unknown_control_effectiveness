function [xfNew,PHIfNew,GfNew,gfNew, PNew, QNew, thetaHatNew,GammaCLNew] = FaultTolUpdateEstimator(x,xf,PHIf,Gf,gf,P,Q,thetaHat,GammaCL,PHI,G,go,scriptYTYStackSum,scriptYxfTStackSum,auxdata,dt)


dxf = (x-xf)/auxdata.kf; %state filter dynamics
dPHIf = (PHI-PHIf)/auxdata.kf; %basis filter dynamics
dGf = (G-Gf)/auxdata.kf; %g*u filter dynamics
dgf = (go-gf)/auxdata.kf; %control effectiveness filter dynamics
clNorm = (1+auxdata.CLNorm*norm(PHIf'*PHIf));
dP = -auxdata.L*P+auxdata.kP1*(PHIf'*PHIf)/clNorm +auxdata.kCL*scriptYTYStackSum; %Least Square P matrix dynamics
dQ = -auxdata.L*Q + auxdata.kP1*PHIf'*((x-xf)/auxdata.kf)/clNorm +auxdata.kCL*scriptYxfTStackSum; %Q matrix  dynamics


M = P*thetaHat-Q; %Error for thetaHat
dthetaHat = -GammaCL*(M); %parameter estimator dynamics
dGammaCL = auxdata.betaCL*GammaCL-GammaCL*P*GammaCL; %Least squares CL gain dynamics

xfNew = xf+dt*dxf; %Update state estimator
PHIfNew = PHIf + dt*dPHIf; %Update filtered basis
GfNew = Gf + dt*dGf; %Update filtered G*u
gfNew = gf + dt*dgf; %Update filtered control effectiveness
PNew = P+dt*dP; %Update Least Squares P matrix
QNew = Q+dt*dQ; %Update Q matrix
thetaHatNew = thetaHat + dt*dthetaHat; %Update parameter estimate
GammaCLNew = GammaCL + dt*dGammaCL; % Update LS CL gain matrix;


%Check thetaHat lies within bounds
[I1] = find(thetaHatNew < auxdata.thetaLower);
[I2] = find(thetaHatNew > auxdata.thetaUpper);

if isempty(I1)~=1;
    for i=1:length(I1)
    thetaHatNew(I1(i)) = auxdata.thetaLower;
    end
end

if isempty(I2)~=1;
    for i=1:length(I2)
    thetaHatNew(I2(i)) = auxdata.thetaUpper;
    end
end

end