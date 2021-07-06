%This file tests out the algorithm to calculate the fault estimates

A = [0,5;0,4;0,0;1,0;1,1]*0;
x = [1;9;3;5;10;20]*2;
B =  [1,2;3,4;6,8]
c = [2;3]
Bx = kron(x,A)

[II,JJ] = find(A(:,:)~=0)
I = II(1);
J = JJ(1);

for i=1:length(x)
    Bx1(i,:)=Bx((i-1)*length(B)+I,J);
    xHat(i,1) = min(max(1/A(I,J)*Bx((i-1)*length(A)+I,J),0),500);
end
Bx1
xHat