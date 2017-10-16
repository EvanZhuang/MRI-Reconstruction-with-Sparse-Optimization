function BPI = myFilteredBackprojectionCentralSlice(sinogram,thetas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered back projection in the frequency domain applying central slice
% theorem-> schlegel & bille 9.2.3
% modified by: Mark Bangert
% m.bangert@dkfz.de 2017

% figure out how big our picture is going to be.
numOfParallelProjections = size(sinogram,1);
numOfAngularProjections  = length(thetas); 

% convert thetas to radians
thetas = (pi/180)*thetas;

% set up the backprojected image
BPI_fourier = complex(zeros(numOfParallelProjections,numOfParallelProjections));

% find the middle index of the projections
midCoord = (numOfParallelProjections+1)/2;

% set up the coords of the image
yCoords = ([1:numOfParallelProjections]) - (numOfParallelProjections+1)/2;

% set up filter
rampFilter = [floor(numOfParallelProjections/2):-1:0 1:ceil(numOfParallelProjections/2-1)];%linspace(-1,1,numOfParallelProjections);

% loop over each projection
for i = 1:numOfAngularProjections

    % figure out which projections to add to which spots
    xCoordsRot = round(midCoord - yCoords*cos(thetas(i)));
    yCoordsRot = round(midCoord - yCoords*sin(thetas(i)));
    
    % convert 2d coords to indices
    ix = sub2ind(numOfParallelProjections*[1 1],xCoordsRot,yCoordsRot); % TRANSPOSE!!!!

    % filter in fourier domain
    filteredProfile_fourier = rampFilter.*fftshift( fft( fftshift(sinogram(:,i)') ) );

    % summation in fourier domain
    BPI_fourier(ix) = BPI_fourier(ix) + filteredProfile_fourier./numOfAngularProjections;

end

% % conversion of aggregated result to spatial domain
% BPI = real( ifftshift(fft2(BPI_fourier)) );
% imagesc(reconstrution2DFT);
% % visualization on the fly
% imagesc(BPI)
BPI = BPI_fourier;
