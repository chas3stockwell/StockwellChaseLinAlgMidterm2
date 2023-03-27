function [ x, niters ] = CG( A, b, x );
        r = b - A * x;
        k = 0; 
        prev_p = 0;

        while ( norm(r) >= eps*norm(b) )
            if (k == 0)
                p = b;
            else
                gamma = (-prev_p'*A*r) / (prev_p'*A*prev_p);
                p = r + gamma * prev_p;
            end
            alpha = p'*r / (p'*A*p);
            x = x + alpha*p;
            r = r - alpha*A*p; 
            prev_p = p; 
            k = k + 1; 
        end
        niters = k;

end