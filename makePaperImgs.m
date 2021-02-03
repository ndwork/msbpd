
function makePaperImgs( mainOut, datacases, testImages )

  outDir = [ mainOut, '/paperImgs' ];

  if ~exist( outDir, 'dir' ), mkdir( outDir ); end
  makeWaveletPlots( outDir );

  nDatacases = numel( datacases );
  for datacaseIndx = 1 : nDatacases
    datacase = datacases( datacaseIndx );
    img = loadDatacase( datacase, testImages );
    sImg = size( img );

    wavSplit = makeWavSplit( sImg );
    
    wtImg = wtDaubechies2( img, wavSplit );
    wtImgScaled = wavScale( abs( wtImg ), wavSplit, 'range', [ 0 1 ] );
    if ~exist( [ outDir, '/wtScaledImg' ], 'dir' ), mkdir( [ outDir, '/wtScaledImg' ] ); end
    imwrite( wtImgScaled, [ outDir, '/wtScaledImg/wtScaledImg_', indx2str(datacase,nDatacases), '.png' ] );

    wavMask = makeLowFreqWavMask( sImg, wavSplit );
    wtImgScaledLow = wtImgScaled( wavMask == 1 );
    wtImgScaledLow = reshape( wtImgScaledLow, [ sqrt(numel(wtImgScaledLow)) sqrt(numel(wtImgScaledLow)) ] );
    wtImgScaledLow = imresize( wtImgScaledLow, 10, 'nearest' );
    if ~exist( [ outDir, '/wtImgScaledLow' ], 'dir' ), mkdir( [ outDir, '/wtImgScaledLow' ] ); end
    imwrite( wtImgScaledLow, [ outDir, '/wtImgScaledLow/wtImgScaledLow_', ...
      indx2str(datacase,nDatacases), '.png' ] );

    acr = makeAutoCalRegion( size(wavMask), wavSplit );
    if ~exist( [ outDir, '/acr' ], 'dir' ), mkdir( [ outDir, '/acr' ] ); end
    imwrite( acr, [ outDir, '/acr/acr_', indx2str(datacase,nDatacases), '.png' ] );
    if ~exist( [ outDir, '/fsr' ], 'dir' ), mkdir( [ outDir, '/fsr' ] ); end
    imwrite( acr>0, [ outDir, '/fsr/fsr_', indx2str(datacase,nDatacases), '.png' ] );

    fftImg = fftshift( ufft2( img ) );
    fftAcrImg = fftImg .* acr;

    lowFreqImg = uifft2( fftAcrImg );
    lowFreqImg = scaleImg( abs(lowFreqImg), [0 1] );
    if ~exist( [ outDir, '/lowFreqImg' ], 'dir' ), mkdir( [ outDir, '/lowFreqImg' ] ); end
    imwrite( lowFreqImg, [ outDir, '/lowFreqImg/lowFreqImg_', indx2str(datacase,nDatacases), '.png' ] );

    fftDiffImg = fftImg - fftAcrImg;
    diffImg = uifft2( ifftshift( fftDiffImg ) );
    wtDiffImg = wtDaubechies2( diffImg, wavSplit );
    wtDiffScaled = wavScale( abs( wtDiffImg ), wavSplit, 'range', [ 0 1 ] );
    if ~exist( [ outDir, '/wtDiffScaled' ], 'dir' ), mkdir( [ outDir, '/wtDiffScaled' ] ); end
    imwrite( wtDiffScaled, [ outDir, '/wtDiffScaled/wtDiffScaled_', indx2str(datacase,nDatacases), '.png' ] );

    wtDiffScaledLow = wtDiffScaled( wavMask == 1 );
    wtDiffScaledLow = reshape( wtDiffScaledLow, [ sqrt(numel(wtDiffScaledLow)) sqrt(numel(wtDiffScaledLow)) ] );
    wtDiffScaledLow = imresize( wtDiffScaledLow, 10, 'nearest' );
    if ~exist( [ outDir, '/wtDiffScaledLow' ], 'dir' ), mkdir( [ outDir, '/wtDiffScaledLow' ] ); end
    imwrite( wtDiffScaledLow, [ outDir, '/wtDiffScaledLow/wtDiffScaledLow_', ...
      indx2str(datacase,nDatacases), '.png' ] );
  end
end


