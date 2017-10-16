function [x, iter, time] = GIST_MCP_Nesterov(A,b,mu,theta,opts)
%solve for min l2-norm + MCP problem
tstart = tic;
fprintf(' ********** start **********\n')
n = size(A,2);
if nargin <= 4
    x0 = zeros(n,1);
    tol = 1e-4;
else
    if isfield(opts,'x0'),      x0      = opts.x0;      else x0 = zeros(n,1);end
    if isfield(opts,'tol'),     tol     = opts.tol;     else tol = 1e-4;     end
end

iter = 0; residue = inf; x = x0; r = residue;
grad = A'*(A*x-b); z = x0; t = 1; beta = max(eig(A*A'));

while (iter <= 1000)
    
    x0 = x;
    y = z - 1/beta * grad;
    %x = sign(y).*max(abs(y)-mu/beta,0);
    x = MCP_Nesterov_subroutine(y, mu, theta);
    t0 = t;
    t = (1 + sqrt(4*t^2 + 1))/2;
    lambda = 1 + (t0 - 1)/t;
    z = x0 + lambda*(x - x0);
    tmp = A*x - b;
    grad = A'*tmp;
    if mod(iter,100) == 0
        r1 = norm(x,1);
        if (r1 <= theta*mu)
            r2 = x'*x/(2*theta);
        else
            r2 = r1 - theta*mu/2;
        end
        r = 0.5*norm(tmp, 2)^2 + mu*(r1 - r2);
        fprintf(' iter = %5d   fval = %12.6f \n', iter, r)
    end
    
    if (norm(x - x0)/max(1,norm(x)) < tol) break; end;
    iter = iter + 1;
end

fprintf(' TERMINATION \n')
time = toc(tstart);
fprintf(' iter = %5d   fval = %12.6f   time = %5f \n', iter+1, r, time)
end
 
