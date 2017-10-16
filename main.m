n0 = 50;

figure
colormap bone 

phantomData = phantom(n0);

phantomData = double(phantomData - min(phantomData(:))); % set min of image to zero

%pad the image with zeros so we don't lose anything when we rotate.
padDims = ceil(norm(size(phantomData)) - size(phantomData));
P       = padarray(phantomData,padDims);

subplot(2,2,1)
imagesc(P);

% set some parameters
freq = 3; % [1/degree]
thetas = 0:1/freq:180-1/freq;

% compute sinogram / radon transformation (schlegel & bille 9.1.1)
sinogram = myRadon(P,thetas);
subplot(2,2,2)
simpleBackprojection = myBackprojection(sinogram,thetas);
reconstrution2DFT = myFilteredBackprojection2DFT(simpleBackprojection);
subplot(2,2,3)
imagesc(reconstrution2DFT);

[c,s]=wavedec2(P,5,'haar');
n = size(c,2); % signal size
m = floor(n/5); % measurement
sigma = 0.05; % noise level

A = randn(m,n);
x0 = c';
e = randn(m,1);
b = A*x0 + sigma*e;
opts.x0 = randn(n,1); opts.tol = 1e-4; mu = 1; theta = 1.2;

for i = 1:10;
   [x, iter, time] = GIST_MCP(A,b,mu,theta,opts);
   mu = mu/2; opts.x0 = x; opts.tol = max(opts.tol/2,1e-6);
end

X_new = waverec2(x',s,'haar');
subplot(2,2,4)
imagesc(X_new)

