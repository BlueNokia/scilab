A = [2 1 -1; -3 -1 2; -2 1 2];
b = [8; -11; -3];

function f = gaussJordan(A,b)
    n = length(b)
    A = [A b]
    for k = 1:n
        for j = n+1:-1:k
            A(k,j) = A(k,j)/A(k,k)
        end
        for i=1:n
            if(i~=k)
                for j = n+1:-1:k
                    A(i,j) = A(i,j)-A(k,j)*A(i,k)
                end
            end
        end
    end
    f = A(:,n+1)
endfunction
disp(gaussJordan(A,b))
