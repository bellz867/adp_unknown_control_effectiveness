function integratedBuffer = FaultToIIntegrateBuffer(timeBuffer,Buffer)
bufferSize = length(timeBuffer);
integratedBuffer = zeros(size(Buffer,1),size(Buffer,2));
for i = 1:bufferSize-1
    integratedBuffer = integratedBuffer + (timeBuffer(i+1)-timeBuffer(i))*(Buffer(:,:,i+1)+Buffer(:,:,i));
end
integratedBuffer = 0.5*integratedBuffer;
end