function makeWaveletPlots( outDir )
  close all

  lineWidth = 3;

%  lpf = [ ...
%    -0.12940952255092145, ...
%    0.22414386804185735, ...
%    0.836516303737469, ...
%    0.48296291314469025, ...
%  ];

lpf = [ ...
  1-sqrt(3), 3-sqrt(3), 3+sqrt(3), 1+sqrt(3) ...
] / ( 4 * sqrt(2) );
%lpf = lpf / norm( lpf );

hpf = [ ...
  -(1+sqrt(3)), 3+sqrt(3), -(3-sqrt(3)), (1-sqrt(3)) ...
] / ( 4 * sqrt(2) );
%hpf = hpf / norm( hpf );

%  hpf = [ ...
%    -0.48296291314469025, ...
%    0.836516303737469, ...
%    -0.22414386804185735, ...
%    -0.12940952255092145, ...
%  ];


  figure;

  subplot(2,1,1);  stemnice( hpf, 'r', 'LineWidth', lineWidth );
  hold on;  stemnice( lpf, 'b', 'LineWidth', lineWidth );
  set( gca, 'XTick', (1:numel(lpf)) );
  title( 'Coefficients' );
  legend( 'High Pass', 'Low Pass' );

  paddedLpf = zeros(512,1);  paddedLpf(1:numel(lpf)) = lpf;
  paddedHpf = zeros(512,1);  paddedHpf(1:numel(hpf)) = hpf;

  fftLpf = 1/sqrt(2) * fftc( paddedLpf );
  fftHpf = 1/sqrt(2) * fftc( paddedHpf );

  fftCoords = size2fftCoordinates( 512 );
  subplot(2,1,2);  plotnice( fftCoords, abs( fftHpf ), 'r', 'LineWidth', lineWidth );
  hold on;  plotnice( fftCoords, abs( fftLpf ), 'b', 'LineWidth', lineWidth );
  title( 'Spectrums' );


  saveas( gcf, [ outDir, '/waveletPlots.png' ] );
  close( gcf );
  
end

