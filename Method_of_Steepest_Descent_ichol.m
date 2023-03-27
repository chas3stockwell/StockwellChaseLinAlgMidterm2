function [x, niters] = Method_of_Steepest_Descent_ichol(A,b,x)
%UNTITLED10 Summary of this function goes here
       L = ichol(sparse(A), struct('type','ict','droptol',1e-3,'michol','off'));
       M = L * L';
       k = 0;
       r = b - A*x; 
       %inputs: specify tolerance. We will use matlab's eps f(x)
       while ( norm(r) >= eps*norm(b) )
            %preconditioning step
            p = M \ r;
            q = A*p; 
            alpha = (p' * r) / (p' * q);
            x = x + alpha * p;
            r = r - alpha * q; 
            k = k + 1;
       end 
       niters = k;


end