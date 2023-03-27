function [ x, niters ]= Method_of_Steepest_Descent( A, b, x )
        %choose a random x. 
       k = 0;
       r = b - A*x; 
       %inputs: specify tolerance. We will use matlab's eps f(x)
       while ( norm(r) >= eps*norm(b) )
            p = r;
            q = A*p; 
            alpha = (p' * r) / (p' * q);
            x = x + alpha * p;
            r = r - alpha * q; 
            k = k + 1;
       end 
       niters = k;
end