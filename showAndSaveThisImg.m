
function showAndSaveThisImg( outDir, img, name, varargin )

  p = inputParser;
  p.addOptional( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.addParameter( 'range', [], @isnumeric );
  p.addParameter( 'showScale', 1, @ispositive );
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.parse( varargin{:} );
  verbose = p.Results.verbose;
  range = p.Results.range;
  showScale = p.Results.showScale;
  wavSplit = p.Results.wavSplit;

  if verbose == true
    figure; imshowscale( img, showScale );
  end

  if numel( wavSplit ) > 0
    outImg = wavScale( img, wavSplit );
  else
    %imshowscale( img, showScale, 'range', range );
    outImg = scaleImg( img, range );
  end

  imwrite( outImg, [ outDir, filesep(), name, '.png' ] );
end
