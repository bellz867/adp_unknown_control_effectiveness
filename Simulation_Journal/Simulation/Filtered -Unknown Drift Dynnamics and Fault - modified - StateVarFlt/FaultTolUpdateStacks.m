function [scriptYTYStackNew,scriptYTxStackNew,minEig,scriptYYTStackSum,scriptYxfTStackSum,stackIndexNew] = FaultTolUpdateStacks(scriptYTYStack,scriptYTxStack,scriptY,scriptG,deltax,thetaHat,stackIndex,auxdata)
% if the last element of the stack has rank higher than 0 then replace every 
% element of stack with the new term to see if it increases min eigenvalue

scriptYTYStackNew = scriptYTYStack;
scriptYTxStackNew = scriptYTxStack;
N = size(scriptYTYStackNew,3); %stack size
stackNorm = (1+auxdata.CLNorm*norm(scriptY'*scriptY));
newYTY = scriptY'*scriptY/stackNorm;
newYTxf = scriptY'*(deltax)/stackNorm;

if (stackIndex <= N) %if the index didnt make it through all the way then the stack is not filled and add the new values at that index
    scriptYTYStackNew(:,:,stackIndex) = newYTY;
    scriptYTxStackNew(:,:,stackIndex) = newYTxf;
    scriptYYTStackSum = sum(scriptYTYStackNew,3); %new Y'Y stack sum
    scriptYxfTStackSum = sum(scriptYTxStackNew,3); %new Y'[x(t)-x(t-T)-Y*thetaHat)] stack sum
    minEig = min(eig(scriptYYTStackSum));%min eigenvalue

else %if it is filled then do the replacements and choose the best combination
    minEigs = zeros(1,N); %min eig for each replacement
    minEigInitial = min(eig(sum(scriptYTYStackNew,3))); %initial min eigenvalue
    scriptYTYStackSumInitial = sum(scriptYTYStackNew,3) + newYTY; %initial sum

    for j = 1:N %do the switch on all stack elements
        scriptYTYSums = scriptYTYStackSumInitial - scriptYTYStackNew(:,:,j); %subtract out the element of the sum and replace with the new one
        minEigs(j) = min(eig(scriptYTYSums)); % find the minimum eigenvalue of the new stack
    end

    [maxMinEig,maxMinEigIndex] = max(minEigs); %get the maximum minimum eigenvalue and its index

    if maxMinEig > minEigInitial %if the max of the replacements has a larger minimum eigenvalue than the initial minimum then take that stack as the new stack
        scriptYTYStackNew(:,:,maxMinEigIndex) = newYTY; %update the Y'Y stack
        scriptYTxStackNew(:,:,maxMinEigIndex) = newYTxf; %update the Y'[x(t)-x(t-T)-Y*thetaHat)] stack
        minEig = maxMinEig;%update the minimimum eigenvalue for the ouptu
    else
        minEig = minEigInitial;
    end
    scriptYYTStackSum = sum(scriptYTYStackNew,3);%set the Y'Y stack sum to the maximizing one
    scriptYxfTStackSum = sum(scriptYTxStackNew,3);%get the sum of the new Y'[x(t)-x(t-T)-Y*thetaHat)] stack
end

stackIndexNew = stackIndex + 1;
end