
function sampleMask = makePoissonDiscSamples4CS( Ns, varargin )

  defaultr = 3 / ( min(Ns) );
  p = inputParser;
  p.addParameter( 'r', defaultr, @isnumeric );
  p.parse( varargin{:} );
  r = p.Results.r;

  kPts = poissonDisc2( r, 'incSize', 25000, 'verbose', true );
  kPts = kPts';

  dks = 1 ./ (Ns-1);
  [~,samples] = movePointsToGrid( kPts, [-0.5, -0.5], 0.5-dks, Ns );
  sampleMask = samples > 0;
  %figure; imshowscale( sampleMask, 3 );

end

