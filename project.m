%  Take an image and project in into a matrix
%  suitable for backprojection
%
%  Calling: ct_data = project (image_data, K)
%
%  Input:  image_data - The gray scale image
%          K          - Number of projections
%
%  Output: ct_data    - The projected data
%
%  version 1.2, 6/11-97 Joergen Arendt Jensen

function ct_data = project (image_data,K)

%  Pre-allocate storage for the projection
%  and ensure that the image is retangular

[N,M]=size(image_data);
N=max(N,M);
image_data(N,N)=0;
ct_data=zeros(N,K);

%  Make the index matrix

xm = ((1:N)'-N/2) * ones(1,M);
ym = ones(N,1) * ((1:M)-M/2);

%  Do the different summations 

dtheta=pi/K;
for i=1:K

  disp([num2str(i/K*100),'% done'])

  %  Rotate the array

  theta = (i-1)*dtheta;
  x = floor( xm*cos(theta) - ym*sin(theta) + N/2 + 0.5);
  y = floor( xm*sin(theta) + ym*cos(theta) + M/2 + 0.5);

  inside = (x>=1) & (x<=N) & (y>=1) & (y<=M);
  x = x.*inside + 1*(1-inside);
  y = y.*inside + 1*(1-inside);

  %  Find the projection

  pro = sum (image_data(x + (y-1)*N) .* inside)';
  ct_data(:,i) = pro(N:-1:1);
  end