
function run_csBetterSampling( varargin )
  close all; rng(1);

  showScale = 1;
  debug = true;
  nIter = 100;
  lambdas = [ 1d-1 1d-2 1d-3 1d-4 5d-2 5d-3 5d-4 5d-5  ];
  vdSigs = [ 100 75 125 150 175 ];
  verbose = false;
  nDatacases = 21;
  noiseSDevs = [ 0 0.008 0.032 0.128 ];
  nSamplesArray = [ 20000 40000 30000 50000 60000 70000 10000 ];
  printEvery = 19;
  mainOut = './out/';
  logFilename = 'log.csv';

  p = inputParser;
  p.addOptional( 'datacases', 0:nDatacases-1, @isnumeric );
  p.parse( varargin{:} );
  datacases = p.Results.datacases;
%datacases = [ 1 2 3 4 8 9 13 7 10 0 11 12 14 15 16 17 18 19 20 ];
datacases = [ 1 2 3 4 8 9 13 7 10 0 ];

  wavSplit = zeros(8);  wavSplit(1,1) = 1;

  testImages = listTestImages();

  if ~exist( mainOut, 'dir' ), mkdir( mainOut ); end
  makePaperImgs( mainOut, wavSplit, datacases, testImages );

  logFile = [ mainOut, filesep, logFilename ];
  if ~exist( logFile, 'file' )
    logID = fopen( logFile, 'a' );
    fprintf( logID, '----  New Run Here ----\n' );
    fprintf( logID, 'noiseSDev, vdSig, nSamples, datacase, Algorithm, lambda, err \n' );
    fclose( logID );
  end

  for noiseSDev = noiseSDevs
    for vdSig = vdSigs
      for nSamples = nSamplesArray

        parfor datacaseIndx = 1 : numel( datacases )
          datacase = datacases( datacaseIndx );
          absDiffRange = [];

          for lambda = lambdas
            lambda_2level = lambda;
            lambda_Nick = lambda;

            outDir = [ mainOut, ...
              '/noiseSDev_', num2str( noiseSDev ), ...
              '/vdSig_', indx2str( vdSig, max( vdSigs ) ), ...
              '/nSamples_', indx2str( nSamples, max( nSamplesArray ) ), ...
              '/datacase_', indx2str( datacase, nDatacases ), ...
              '/lambda_', num2str( lambda ), '/' ];
            if ~exist( outDir, 'dir' ), mkdir( outDir ); end

            %[img,lambda,lambda_2level,lambda_Nick] = loadDatacase( datacase );
            img = loadDatacase( datacase, testImages );
            sImg = size( img );

            disp( [ 'Working on ', outDir ] );
            if exist( [ outDir, filesep, 'absDiff_2level.png' ], 'file' )
              disp( '  Previously completed.  Continuing.' );
              continue;
            end

            showAndSaveThisImg( outDir, abs(img), 'origImg', verbose, 'showScale', showScale );

            wtImg = wtDeaubechies2( img, wavSplit );
            showAndSaveThisImg( outDir, abs(wtImg), 'absWtImg', verbose, 'showScale', showScale, ...
              'wavSplit', wavSplit );
            fftImg = fftshift( ufft2( img ) );

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

            for i = 1:nRecons
            %parfor i = 1:nRecons
              if i == 1
                [ thisRecon, theseOValues ] = csReconFISTA( fftSamples, lambda, 'wavSplit', wavSplit, ...
                  'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 2
                [ thisRecon, theseOValues ] = csReconFISTA_nick( fftSamples, lambda_Nick, 'wavSplit', wavSplit, ...
                  'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 3
                [ thisRecon, theseOValues ] = csReconFISTA( fftSamples_wACR, lambda, 'wavSplit', wavSplit, ...
                  'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 4
                [ thisRecon, theseOValues ] = csReconFISTA_nick( fftSamples_wACR, lambda_Nick, ...
                  'wavSplit', wavSplit, 'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );

              elseif i == 5
                [ thisRecon, theseOValues, xStar ] = csReconFISTA_2level( fftSamples_wACR, lambda_2level, ...
                  'wavSplit', wavSplit, 'verbose', true, 'printEvery', printEvery, 'debug', debug, 'nIter', nIter );
                showAndSaveThisImg( outDir, abs(xStar), 'absXStar', verbose, 'showScale', showScale, ...
                  'wavSplit', wavSplit, 'range', abs(wtImg) );

              end

              thisRecon = abs( thisRecon );

              recons{i} = thisRecon;
              oValues{i} = theseOValues;
            end
            recon = recons{1};
            recon_nick = recons{2};
            recon_wACR = recons{3};
            recon_nick_wACR = recons{4};
            recon_2level = recons{5};

            oValues_nick = oValues{2};
            oValues_wACR = oValues{3};
            oValues_nick_wACR = oValues{4};
            oValues_2level = oValues{5};
            oValues = oValues{1};

            showAndSaveThisPlot( outDir, oValues, 'oValues', verbose );
            showAndSaveThisPlot( outDir, oValues_nick, 'oValues_nick', verbose );
            showAndSaveThisPlot( outDir, oValues_wACR, 'oValues_wACR', verbose );
            showAndSaveThisPlot( outDir, oValues_nick_wACR, 'oValues_nick_wACR', verbose );
            showAndSaveThisPlot( outDir, oValues_2level, 'oValues_2level', verbose );

            showAndSaveThisImg( outDir, abs(img), 'origImg', verbose, 'showScale', showScale );
            showAndSaveThisImg( outDir, abs(recon), 'csRecon', verbose, 'showScale', showScale, ...
              'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_nick), 'csReconNick', verbose, 'showScale', showScale, ...
              'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_wACR), 'csRecon_wACR', verbose, 'showScale', showScale, ...
              'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_nick_wACR), 'csRecon_nick_wACR', verbose, 'showScale', showScale, ...
              'range', abs(img) );
            showAndSaveThisImg( outDir, abs(recon_2level), 'csRecon2Level', verbose, 'showScale', showScale, ...
              'range', abs(img) );

            logID = fopen( [ mainOut, filesep, logFilename ], 'a' );

            diff = recon - img;
            err = norm( diff(:) ) / norm( img(:) );
            disp([ 'Err: ', num2str( err ) ]);

            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err, ', ...
              num2str( lambda ), ', ', ...
              num2str(err), ...
              '\n' ] );

            diff_nick = recon_nick - img;
            err_nick = norm( diff_nick(:) ) / norm( img(:) );
            disp([ 'Err Nick: ', num2str( err_nick ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_nick, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_nick ), ...
              '\n' ] );

            diff_wACR = recon_wACR - img;
            err_wACR = norm( diff_wACR(:) ) / norm( img(:) );
            disp([ 'Err wACR: ', num2str( err_wACR ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_wACR, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_wACR ), ...
              '\n' ] );

            diff_nick_wACR = recon_nick_wACR - img;
            err_nick_wACR = norm( diff_nick_wACR(:) ) / norm( img(:) );
            disp([ 'Err Nick_wACR: ', num2str( err_nick_wACR ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_nick_wACR, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_nick_wACR ), ...
              '\n' ] );

            diff_2level = recon_2level - img;
            err_2level = norm( diff_2level(:) ) / norm( img(:) );
            disp([ 'Err 2level: ', num2str( err_2level ) ]);
            fprintf( logID, [ ...
              num2str( noiseSDev ), ', ', ...
              indx2str( vdSig, max( vdSigs ) ), ', ', ...
              indx2str( nSamples, max( nSamplesArray ) ), ', ', ...
              indx2str( datacase, nDatacases ), ', ', ...
              'Err_2level, ', ...
              num2str( lambda ), ', ', ...
              num2str( err_2level ), ...
              '\n' ] );

            if numel( absDiffRange ) == 0
              absDiffRange = [ 0 max( abs(diff(:) ) ) ];
            end

            showAndSaveThisImg( outDir, abs(diff), 'absDiff', verbose, 'showScale', showScale );
            showAndSaveThisImg( outDir, abs(diff_wACR), 'absDiff_wACR', verbose, 'showScale', showScale, ...
              'range', absDiffRange );
            showAndSaveThisImg( outDir, abs(diff_nick), 'absDiff_nick', verbose, 'showScale', showScale, ...
              'range', absDiffRange );
            showAndSaveThisImg( outDir, abs(diff_nick_wACR), 'absDiff_nick_wACR', verbose, 'showScale', showScale, ...
              'range', absDiffRange );
            showAndSaveThisImg( outDir, abs(diff_2level), 'absDiff_2level', verbose, 'showScale', showScale, ...
              'range', absDiffRange );

            fclose( logID );
          end
        end
      end
    end
  end

end
