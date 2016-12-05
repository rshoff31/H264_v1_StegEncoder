function z = mirror_pad(z)

[n,m] = size(z);

for i = 1:n
    for j = 2:m-1
        
        if z(i,j+1) == 0 & z(i,j) ~= 0
            z(i,j+1) = z(i,j);
        end
    end
    for j = m:-1:2

        if z(i,j-1) == 0 & z(i,j) ~= 0
            z(i,j-1) = z(i,j);
        end
    end
end
for j = 1:m
    for i = 2:n-1
        
        if z(i+1,j) == 0 & z(i,j) ~= 0
            z(i+1,j) = z(i,j);
        end
    end
    for i = n:-1:2
        
        if z(i-1,j) == 0 & z(i,j) ~= 0
            z(i-1,j) = z(i,j);
        end
    end
end
