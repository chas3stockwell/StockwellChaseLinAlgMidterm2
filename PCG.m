function [x,niters] = PCGa(A,b,x)
   L = ichol(sparse(A), struct('type','ict','droptol',1e-3,'michol','off'));
   M = L * L';
   k = 0;
   r = b - A * x;
   
   prev_r = zeros(size(A, 2), 1);
   prev_p = zeros(size(A, 2), 1);
   
   while ( norm(r) >= eps*norm(b) )
        z = M \ r;
        if k == 0
            p = z;
        else    
            gamma = (r'* (z) ) / (prev_r' * (M \ prev_r) );
            %gamma = (z' * r) / ((M \ prev_r)' * prev_r);
            p = z + gamma*prev_p;
        end
        q = A*p;
        alpha = (z'*r) / (p'*q);
        
        prev_p = p; 
        prev_r = r;
        x = x + alpha*p;
        r = r - alpha*q;

        k = k + 1; 
    end
    niters = k;  
end
