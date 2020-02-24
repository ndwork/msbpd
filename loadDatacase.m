
function [img,lambda,lambda_2level,lambda_Nick] = loadDatacase( datacase, testImages )

  lambda = 1d-3;
  lambda_Nick = lambda;
  lambda_2level = lambda;

  if datacase == 0
    load( './mrData/brainData.mat', 'brainData' );
    padded = padData( brainData, [512 512] );
    img = ifftc( padded );
    img = scaleImg( img, [0 1], [0 4000] );
  end

  if datacase == 1
    load( './mrData/mriImg.mat', 'mriImg' );
    img = scaleImg( mriImg, [0 1], [0 0.7] );

  elseif datacase == 2
    load( './mrData/nian1.mat', 'nian1' );
    img = scaleImg( nian1, [0 1] );
    
  elseif datacase == 3
    load( './mrData/nian2.mat', 'nian2' );
    img = scaleImg( nian2, [0 1] );

  elseif datacase > 3
    testImgIndx = datacase - 3;
    testImgName = testImages( testImgIndx ).name;
    img = double( imread( testImgName ) ) / 255.;
    if ~ismatrix( img ), img = rgb2gray( img ); end

  end
  
  img = imresize( img, [512 512] );
end