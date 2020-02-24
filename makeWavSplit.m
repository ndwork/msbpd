
function wavSplit = makeWavSplit( sImg )

  binPow = logBase( sImg, 2 );

  if max( mod( binPow(:), 1 ) ) ~= 0
    error( 'Image must have sizes of power of 2' );
  end

  wavSplit = zeros( sImg ./ 64 );
  wavSplit(1) = 1;
end
