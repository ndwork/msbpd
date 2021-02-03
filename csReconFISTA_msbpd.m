

function [recon,oValues,xStar] = csReconFISTA_msbpd( samples, lambda, varargin )
  % recon = csReconFISTA_msbpd( samples, lambda [, 'debug', debug, ...
  %   'nIter', nIter, 'printEvery', printEvery, 'wavSplit', wavSplit, ...
  %   'verbose', verbose, 'waveletType', waveletType ] )
  %
  % This routine minimizes 0.5 * || Ax - b ||_2^2 + lambda || W x ||_1
  %   where A is sampleMask * Fourier Transform * real part, and
  %   W is the wavelet transform.  It assumes a fully sampled center region
  %   corresponding the the two-level sampling scheme defined in
  %   "Breaking the coherence barrier: A new theory for compressed sensing"
  %   by Adcock, Ben, et al.
  %   Moreover, it uses the fully sampled scheme in a way to increase
  %   sparsity in the optimiation problem
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % debug - if true, reduces the default number of iterations to 30 and forces verbose
  %         statements during optimization
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % polish - if set to true, adds a polishing step (default is false)
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % verbose - if true, prints informative statements
  % waveletType - either 'Daubechies' for Daubechies-4 (default) or 'Haar'
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  wavSplit = zeros(8);  wavSplit(1,1) = 1;

  p = inputParser;
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'debug', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'polish', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'waveletType', 'Daubechies', @(x) true );
  p.parse( varargin{:} );
  checkAdjoints = p.Results.checkAdjoints;
  debug = p.Results.debug;
  nIter = p.Results.nIter;
  printEvery = p.Results.printEvery;
  wavSplit = p.Results.wavSplit;
  waveletType = p.Results.waveletType;
  verbose = p.Results.verbose;

  if numel( nIter ) == 0
    if debug == true
      nIter = 30;
    else
      nIter = 100;
    end
  end

  M = ( samples ~= 0 );
  sCenterRegion = size( samples ) ./ ( 2 .^ ( log2( size( wavSplit ) ) + 1 ) );
  samplesL = zeroOuterRegion( samples, sCenterRegion );
  
  acr = makeAutoCalRegion( size( samples ), wavSplit );
  
  acrSamplesL = acr .* samplesL;
  beta = samples( M == 1 ) - acrSamplesL( M == 1 );

  % RI = [ Re; Im; ]
  % A = M F RI , A' = (RI)' * F' * M
  % A' * A = (RI)' * F' * M * F * RI
  % gGrad = A'*A*x - A'*b;

  if strcmp( waveletType, 'Daubechies' )
    W = @(x) wtDaubechies2( x, wavSplit );
    WT = @(y) iwtDaubechies2( y, wavSplit );
  elseif strcmp( waveletType, 'Haar' )
    W = @(x) wtHaar2( x, wavSplit );
    WT = @(y) iwtHaar2( y, wavSplit );
  else
    error( 'Unrecognized wavelet type' );
  end

  function out = F( x )
    out = fftshift( fftshift( ufft2( ifftshift( ifftshift( x, 1 ), 2 ) ), 1 ), 2 );
  end

  function out = Fadj( y )
    out = fftshift( fftshift( uifft2( ifftshift( ifftshift( y, 1 ), 2 ) ), 1 ), 2 );
  end

  function out = A( x )
    WTx = WT( x );
    FWTx = F( WTx );
    out = FWTx( M == 1 );
  end

  function out = Aadj( x )
    MTx = zeros( size(M) );
    MTx( M == 1 ) = x;
    FadjMTx = Fadj( MTx );
    out = W( FadjMTx );
  end

  if checkAdjoints == true
    % Variable used during debugging of this routine
    x1 = rand( size(samples) ) + 1i * rand( size( samples ) );
    if checkAdjoint( x1, W, WT ) ~= true, error( 'WT is not the transpose of W' ); end
    if checkAdjoint( x1, @F, @Fadj ) ~= true, error( 'Fadj is not the transpose of F' ); end
    if checkAdjoint( x1, @A, @Aadj ) ~= true, error( 'Aadj is not the transpose of A' ); end
  end

  function out = g( x )
    diff = A( x ) - beta;
    out = 0.5 * norm( diff(:), 2 ).^2;
  end

  AadjBeta = Aadj( beta );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - AadjBeta;
  end

  nPixels = numel( samples );
  proxth = @(x,t) proxL1Complex( x, t * lambda / nPixels );

  function out = h( x )
    out = sum( abs( x(:) ) ) * lambda / nPixels;
  end

  x0 = W( Fadj( samples ) );
  wavMask = makeLowFreqWavMask( size( samples ), wavSplit );
  x0 = x0 .* (1-wavMask);

  t = 1;
  if debug
    %[xStar,oValues] = fista( x0, @g, @gGrad, proxth, 'h', @h, 'verbose', verbose );   %#ok<ASGLU>
    [xStar,oValues] = fista_wLS( x0, @g, @gGrad, proxth, 'h', @h, ...
      't0', t, 'N', nIter, 'verbose', true, 'printEvery', printEvery );                                                                     %#ok<ASGLU>
  else
    %xStar = fista( x0, @g, @gGrad, proxth );   %#ok<UNRCH>
    xStar = fista_wLS( x0, @g, @gGrad, proxth, 't0', t, 'N', nIter, ...
      'verbose', verbose, 'printEvery', printEvery );
    oValues = [];
  end

  %recon = WT( xStar ) + Fadj( samplesL );
  recon = WT( xStar ) + Fadj( acrSamplesL );
end

