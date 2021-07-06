function v = SimplexVert(n)
% This function calculates the (n+1) vertices of an n-dimensional simplex 


%Ouput is (n) x (n+1)
c = zeros(n,n+1);

for i = 1:n
    c(i,i) = sqrt(1-sum(c(1:i-1,i).^2));
    
    for j = i+1:n+1
        c(i,j) = (-1/n-(c(1:i-1,i)'*c(1:i-1,j)))/c(i,i);
       
    end
end

v = flipud(c);

end