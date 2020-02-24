
function showAndSaveThisPlot( outDir, oValues, name, varargin )

  p = inputParser;
  p.addOptional( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  verbose = p.Results.verbose;

  figH = figure;
  semilogynice( oValues );
  titlenice( name );
  saveas( figH, [ outDir, filesep(), name, '.png' ] );
  if verbose ~= true, close( figH ); end

end
