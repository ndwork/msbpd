
function out = uifftc( in )
  out = fftc( in ) * sqrt( numel( in ) );
end
