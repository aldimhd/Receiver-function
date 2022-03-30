% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = gfilter( x, nfft, gauss, dt)
%
%    convolve a function with a unit-area Gaussian filter.
%  
% Liguria & Ammon, use G(omega)=exp(-omega^2/(4.*L^2))  
%  G(f0)=exp(-1/2)=.606
% 
%

% get signal in fourier domain
Xf = fft( x, nfft ); 
%Xf(0) is real, Real(Xf(1)) = Real(Xf(nfft))

% Convolve signal with filter, the dt makes the same units
Xf = Xf.*gauss*dt;

% Back to time domain
xnew = real( ifft( Xf, nfft) );