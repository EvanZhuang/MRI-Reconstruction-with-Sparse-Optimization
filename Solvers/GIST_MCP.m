function [x, iter, time] = GIST_MCP(A,b,mu,theta,opts)
%solve for min l2-norm + MCP problem
tstart = tic; 
fprintf(' ********** start **********\n')
n = size(A,2); sigma = 1e-4;
M = 5; % nonmonotone parameter
%g = @(x) At*A*x-At*b;
if nargin <= 4
  x0 = zeros(n,1);
  tol = 1e-4;
else
  if isfield(opts,'x0'),      x0      = opts.x0;      else x0 = zeros(n,1);end
  if isfield(opts,'tol'),     tol     = opts.tol;     else tol = 1e-4;     end
end
iter = 1; x = x0;
r1 = norm(x,1); 
I = [r1 <= theta*mu];
tmp = A*x - b;
fval = 0.5*norm(tmp, 2)^2 + mu*norm(x(I),1) - norm(x(I))^2/2/theta + (sum(1 - I))*theta*mu^2/2;
grad = A'*tmp;

f1 = -inf*ones(M,1);
f1(1) = fval;
fmax = max(f1);


while 1 == 1
    if (iter > 1)
        w = x - x0;  y = grad - grad0;
        bbstep =max(min((w'*y)/(w'*w),1e8),1e-8);
    else
        bbstep = 1;
    end
    x0 = x; grad0 = grad;
    x = MCPsubroutine(x0-1/bbstep*grad0, mu, theta); % change input to (x-1/bbstep*grad, mu, theta)
    r1 = norm(x,1);
    I = [r1 <= theta*mu];
    tmp = A*x - b;
    fnew = 0.5*norm(tmp, 2)^2 + mu*norm(x(I),1) - norm(x(I))^2/2/theta + (sum(1 - I))*theta*mu^2/2;
    while (fnew >= fmax - sigma/2*norm(x - x0)^2)
        x = MCPsubroutine(x0-1/bbstep*grad0, mu, theta);
        tmp = A*x - b;
        bbstep = bbstep*2;
        r1 = norm(x,1);
        I = [r1 <= theta*mu];
        fnew = 0.5*norm(tmp, 2)^2 + mu*norm(x(I),1) - norm(x(I))^2/2/theta + (sum(1 - I))*theta*mu^2/2;
    end
    grad = A'*tmp;
    f1(mod(iter,M)+1) = fnew;
    fmax = max(f1);
    if (norm(x - x0)/max(1,norm(x)) < tol)
        break
    end
    iter = iter + 1;
    if (mod(iter,100) == 0); fprintf(' iter = %5d   fval = %12.6f \n', iter, fnew)
    end
end

fprintf(' TERMINATION \n')
time = toc(tstart);
fprintf(' iter = %5d   fval = %12.6f   time = %5f \n', iter, fnew, time)
 
