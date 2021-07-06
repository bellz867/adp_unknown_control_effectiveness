function muHatNew = FaultTolGetFaultEstimate(WmuHat,muHatOld,m,auxdata)


n = size(WmuHat,2);

muHatNew = zeros(m,n); %Initialize the fault estimate   
    for r=1:n
        I = r;  %Get the rth entry for the row
        J = r;  %Get the rth entry for the column        
        for k=1:m
            muHatNew(k,r) =-WmuHat((k-1)*n+I,J); %Update if WgHat is non-zero
        end
    end
        muHatNew = sum(muHatNew,2)/n; %Get Average across each column

end
    