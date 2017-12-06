% This tests the idea of compressed sensing

m = 300; % measurement
n = 5000; % signal size
s = 50;
sigma = 0.05; % noise level

A = randn(m,n);
x0 = zeros(n,1);
I = randperm(n);
x0(I(1:s)) = randn(s,1);
e = randn(m,1);
b = A*x0 + sigma*e;
opts.x0 = zeros(n,1); opts.tol = 1e-4; mu = 1; theta = 1.2;

profile on
for i = 1:10;
   [x, iter, time] = GIST_MCP_Nesterov(A,b,mu,theta,opts);
   mu = mu/2; opts.x0 = x; opts.tol = max(opts.tol/2,1e-6);
end
profile off

plot(x0,'ro');
hold on
plot(x,'b*');
