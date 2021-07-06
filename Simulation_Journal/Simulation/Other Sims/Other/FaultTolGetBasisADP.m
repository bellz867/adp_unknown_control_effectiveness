function [sig, sigPrime] = FaultTolGetBasisADP(y, auxdata)
v2 = auxdata.scale;
db = auxdata.dbar;
D = auxdata.CenterOffsetsScale*auxdata.CenterOffsets;
n = length(y);
sig = zeros(auxdata.nodes,1);
sigPrime = zeros(auxdata.nodes,n);

x = y./auxdata.NormDyn; %Normalize states

if (auxdata.BasisType == 1)
    for i=1:size(D,2)
        if (auxdata.ShrinkingCenterOffset == 1)
            nm = (x'*x+db)./(1+v2*(x'*x));
            sig(i,1) = x'*(x+nm*D(:,i));
            sigPrime(i,:) = (2*x+nm.*D(:,i)+2*(1-v2*db)*x*D(:,i)'*x/(1+v2*(x'*x))^2)';
        else 
            nm = auxdata.CenterOffsetsScale;
            sig(i,1) = x'*(x+nm*D(:,i));
            sigPrime(i,:) = (2*x+nm.*D(:,i))';
        end
        
    end
    
else
    for i=1:size(D,2)
        if (auxdata.ShrinkingCenterOffset == 1)
            nm = (x'*x+db)./(1+v2*(x'*x));
            sig(i,1) = exp(x.'*(x+nm.*D(:,i)))-1;
            sigPrime(i,:)=(2*x+nm.*D(:,i)+2*(1-v2*db)*x*D(:,i)'*x/(1+v2*(x.'*x))^2)'*exp(x'*(x+nm.*D(:,i)));
        else
            nm = auxdata.CenterOffsetsScale;
            sig(i,1) = exp(x.'*(x+nm.*D(:,i)))-1;
            sigPrime(i,:)=(2*x+nm.*D(:,i))'*exp(x'*(x+nm.*D(:,i)));
        end
    end
    





end