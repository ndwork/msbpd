
function run_msbpd( varargin )
  close all; rng(1);

  showScale = 1;
  debug = true;
  nIter = 100;
  lambdas = [ 1d5 1d4 1d3 1d2 1d1 1d0 1d-1 1d-2 1d-3 1d-4 5d1 5d2 5d3 5d4 5d5 2d1 2d2 2d3 2d4 2d5 ];
  vdSigs = [ 100 75 125 150 ];
  verbose = false;
  nDatacases = 21;
  noiseSDevs = [ 0 0.008 0.032 0.128 ];
  nSamplesArray = [ 20000 40000 30000 50000 60000 70000 10000 ];
  printEvery = 19;
  mainOut = './out/';
  logFilename = 'log.csv';

  datacases = [ 1 4 8 9 13 ];
  %datacases = 0 : nDatacases-1;
  p = inputParser;
  p.addOptional( 'datacases', datacases, @isnumeric );
  p.parse( varargin{:} );
  datacases = p.Results.datacases;


  testImages = listTestImages();

  if ~exist( mainOut, 'dir' ), mkdir( mainOut ); end
  makePaperImgs( mainOut, datacases, testImages );

  for noiseSDev = noiseSDevs
    for vdSig = vdSigs
      for nSamples = nSamplesArray

        %parfor datacaseIndx = 1 : numel( datacases )
        for datacaseIndx = 1 : numel( datacases )
          datacase = datacases( datacaseIndx );
          absDiffRange = [];

          thisOut = [ mainOut, filesep, 'noiseSDev_', num2str( noiseSDev ), ...
            filesep, 'vdSig_', indx2str( vdSig, max( vdSigs ) ), ...
            filesep, 'nSamples_', indx2str( nSamples, max( nSamplesArray ) ), ...
            filesep, 'datacase_', indx2str( datacase, nDatacases ) ];
          if ~exist( thisOut, 'dir' ), mkdir( thisOut ); end

          logFile = [ thisOut, filesep, logFilename ];
          if ~exist( logFile, 'file' )
            logID = fopen( logFile, 'a' );
            fprintf( logID, 'noiseSDev, vdSig, nSamples, datacase, Algorithm, lambda, err, ssim \n' );
            fclose( logID );
          end
          
          for lambda = lambdas
            lambda_msbpd = lambda;
            lambda_maskLF = lambda;

            outDir = [ thisOut, ...
              filesep, 'lambda_', num2str( lambda ), filesep ];
            if ~exist( outDir, 'dir' ), mkdir( outDir ); end

            %[img,lambda,lambda_msbpd,lambda_maskLF] = loadDatacase( datacase );
            img = loadDatacase( datacase, testImages );
            sImg = size( img );

            wavSplit = makeWavSplit( sImg );
            
            disp( [ 'Working on ', outDir ] );
            if exist( [ outDir, filesep, 'absDiff_msbpd.png' ], 'file' )
              disp( '  Previously completed.  Continuing.' );
              continue;
            end

            showAndSaveThisImg( outDir, abs(img), 'origImg', verbose, 'showScale', showScale );

            wtImg = wtDaubechies2( img, wavSplit );
            showAndSaveThisImg( outDir, abs(wtImg), 'absWtImg', verbose, 'showScale', showScale, ...
              'wavSplit', wavSplit );
            fftImg = fftshift( fftshift( ufft2( ifftshift( ifftshift( img, 1 ), 2 ) ), 1 ), 2 );

            acr = makeAutoCalRegion( sImg, wavSplit );

            %r = 1.25 / ( min(sImg) );
            %sampleMask = makePoissonDiscSamples4CS( sImg, 'r', r );
            sampleMask = vdSampleMask( sImg, vdSig, nSamples );
            sampleMask_wACR = sampleMask | acr;
            showAndSaveThisImg( outDir, sampleMask_wACR, 'sampleMask_wACR', verbose, ...
              'showScale', showScale );
            nSamples_wACR = sum( sampleMask_wACR(:) );
            fftSamples_wACR = sampleMask_wACR .* fftImg;
            fftSamples_wACR = fftSamples_wACR + noiseSDev * ( randn( sImg ) + 1i * randn( sImg ) );

            searchF = @(in) sum( sum( vdSampleMask( sImg, vdSig, round( in ) ) ) ) - nSamples_wACR;
            nSamples2Make = binarySearch( searchF, 1, nSamples_wACR*10, 'tol', 0.5 );
            sampleMask = vdSampleMask( sImg, vdSig, round( nSamples2Make ) );
            nSamplesInMask = sum( sampleMask(:) );
            showAndSaveThisImg( outDir, sampleMask, 'sampleMask', verbose, 'showScale', showScale );
            fftSamples = sampleMask .* fftImg;
            fftSamples = fftSamples + noiseSDev * ( randn( sImg ) + 1i * randn( sImg ) );

            disp([ 'nSamples: ', num2str(nSamplesInMask), '   nSamples_wACR: ', num2str(nSamples_wACR) ]);

            nRecons = 5;
            recons = cell( nRecons, 1 );
            oValues = cell( nRecons, 1 );

            parfor i = 1:nRecons
              if i == 1
                [ thisRecon, theseOValues ] = csReconFISTA( fftSamples, lambda, 'wavSplit', wavSplit, ...
                  'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 2
                [ thisRecon, theseOValues ] = csReconFISTA_maskLF( fftSamples, lambda_maskLF, 'wavSplit', wavSplit, ...
                  'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 3
                [ thisRecon, theseOValues ] = csReconFISTA( fftSamples_wACR, lambda, 'wavSplit', wavSplit, ...
                  'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 4
                [ thisRecon, theseOValues ] = csReconFISTA_maskLF( fftSamples_wACR, lambda_maskLF, ...
                  'wavSplit', wavSplit, 'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 5
                [ thisRecon, theseOValues, xStar ] = csReconFISTA_msbpd( fftSamples_wACR, lambda_msbpd, ...
                  'wavSplit', wavSplit, 'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );
                showAndSaveThisImg( outDir, abs(xStar), 'absXStar', verbose, 'showScale', showScale, ...
                  'wavSplit', wavSplit, 'range', abs(wtImg) );

              end

              thisRecon = abs( thisRecon );

              recons{i} = thisRecon;
              oValues{i} = theseOValues;
            end
            recon = recons{1};
            recon_maskLF = recons{2};
            recon_wACR = recons{3};
            recon_maskLF_wACR = recons{4};
            recon_msbpd = recons{5};

            oValues_maskLF = oValues{2};
            oValues_wACR = oValues{3};
            oValues_maskLF_wACR = oValues{4};
            oValues_msbpd = oValues{5};
            oValues = oValues{1};

            showAndSaveThisPlot( outDir, oValues, 'oValues', verbose );
            showAndSaveThisPlot( outDir, oValues_maskLF, 'oValues_maskLF', verbose );
            showAndSaveThisPlot( outDir, oValues_wACR, 'oValues_wACR', verbose );
            showAndSaveThisPlot( outDir, oValues_maskLF_wACR, 'oValues_maskLF_wACR', verbose );
            showAndSaveThisPlot( outDir, oValues_msbpd, 'oValues_msbpd', verbose );

            showAndSaveThisImg( outDir, abs(img), 'origImg', verbose, 'showScale', showScale );
            showAndSaveThisImg( outDir, abs(recon), 'csRecon', verbose, 'showScale', showScale, ...
              'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_maskLF), 'csReconFISTA_maskLF', verbose, ...
              'showScale', showScale, 'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_wACR), 'csRecon_wACR', verbose, 'showScale', showScale, ...
              'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_maskLF_wACR), 'csRecon_maskLF_wACR', verbose, ...
              'showScale', showScale, 'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_msbpd), 'csRecon_msbpd', verbose, ...
              'showScale', showScale, 'range', abs(img) );

            logID = fopen( [ thisOut, filesep, logFilename ], 'a' );

            diff = recon - img;
            err = norm( diff(:) ) / norm( img(:) );
            ssimValue = ssim( recon, img );
            disp([ 'Err: ', num2str( err ) ]);
            disp([ 'ssimValue: ', num2str( ssimValue ) ]);

            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err, ', ...
              num2str( lambda ), ', ', ...
              num2str(err), ', ', ...
              num2str(ssimValue), ...
              '\n' ] );

            diff_maskLF = recon_maskLF - img;
            err_maskLF = norm( diff_maskLF(:) ) / norm( img(:) );
            ssimValue_maskLF = ssim( recon_maskLF, img );
            disp([ 'Err Nick: ', num2str( err_maskLF ) ]);
            disp([ 'ssimValue_maskLF: ', num2str( ssimValue_maskLF ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_maskLF, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_maskLF ), ', ', ...
              num2str( ssimValue_maskLF ), ...
              '\n' ] );

            diff_wACR = recon_wACR - img;
            err_wACR = norm( diff_wACR(:) ) / norm( img(:) );
            ssimValue_wACR = ssim( recon_wACR, img );
            disp([ 'Err wACR: ', num2str( err_wACR ) ]);
            disp([ 'ssimValue_wACR: ', num2str( ssimValue_wACR ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_wACR, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_wACR ), ', ', ...
              num2str( ssimValue_wACR ), ...
              '\n' ] );

            diff_maskLF_wACR = recon_maskLF_wACR - img;
            err_maskLF_wACR = norm( diff_maskLF_wACR(:) ) / norm( img(:) );
            ssimValue_maskLF_wACR = ssim( recon_maskLF_wACR, img );
            disp([ 'Err Nick_wACR: ', num2str( err_maskLF_wACR ) ]);
            disp([ 'ssimValue_maskLF_wACR: ', num2str( ssimValue_maskLF_wACR ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_maskLF_wACR, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_maskLF_wACR ), ', ', ...
              num2str( ssimValue_maskLF_wACR ), ...
              '\n' ] );

            diff_msbpd = recon_msbpd - img;
            err_msbpd = norm( diff_msbpd(:) ) / norm( img(:) );
            ssimValue_msbpd = ssim( recon_msbpd, img );
            disp([ 'Err msbpd: ', num2str( err_msbpd ) ]);
            disp([ 'SSIM msbpd: ', num2str( ssimValue_msbpd ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_msbpd, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_msbpd ), ', ', ...
              num2str( ssimValue_msbpd ), ...
              '\n' ] );

            fclose( logID );

            if numel( absDiffRange ) == 0
              absDiffRange = [ 0 max( abs(diff(:) ) ) ];
            end

            showAndSaveThisImg( outDir, abs(diff), 'absDiff', verbose, 'showScale', showScale );
            showAndSaveThisImg( outDir, abs(diff_wACR), 'absDiff_wACR', verbose, 'showScale', showScale, ...
              'range', absDiffRange );
            showAndSaveThisImg( outDir, abs(diff_maskLF), 'absDiff_maskLF', verbose, 'showScale', showScale, ...
              'range', absDiffRange );
            showAndSaveThisImg( outDir, abs(diff_maskLF_wACR), 'absDiff_maskLF_wACR', verbose, ...
              'showScale', showScale, 'range', absDiffRange );
            showAndSaveThisImg( outDir, abs(diff_msbpd), 'absDiff_msbpd', verbose, 'showScale', showScale, ...
              'range', absDiffRange );
          end
        end
      end
    end
  end

end
