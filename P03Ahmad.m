function [b] = P03Ahmad(L,l,a)
[n,m] = size(a);
k = 1;
N(k) = ceil(sum(a)/L) + 1;

v = repmat(zeros(n,1), [N(k)-1 1]);
c = [v;a];

v = repmat(L,[N(k)-1 1]);
b = [l;v];
beq  = ones(n,1);


for p = 1:N(k)
    A(p,(p-1)*n+1:p*n) = a;
end

for q = 1:n
    for r = 1:N(k)
        Aeq(q,q+(r-1)*n) = 1;
    end
end

x = bintprog(c,A,b,Aeq,beq);
min = c'*x;
if (min == 0)
    k = k+1;
    N(k) = N(k-1)-1;
    clear x;
    clear c;
    clear A;
    clear b;
    clear Aeq;
    clear beq;
    v = repmat(zeros(n,1), [N(k)-1 1]);
    c = [v;a];
    
    v = repmat(L,[N(k)-1 1]);
    b = [l;v];
    beq  = ones(n,1);
    
    
    for p = 1:N(k)
        A(p,(p-1)*n+1:p*n) = a;
    end
    
    for q = 1:n
        for r = 1:N(k)
            Aeq(q,q+(r-1)*n) = 1;
        end
    end
    
    
    
    x = bintprog(c,A,b,Aeq,beq);
end





x = reshape(x,n,N(k));

w = zeros(n,N(k));
for j = 1:N(k)
    for i = 1:n
        if(x(i,j) == 1)
            w(i,j) = a(i);
        end
    end
end

w = reshape(w,n*N(k),1);
b = nonzeros(w);
end
       