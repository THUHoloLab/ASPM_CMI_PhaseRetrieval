%ASM prepared for the waveguide imaging
% distance scalar: mm


function I_ASM=ASM_diffraction(lambda,I0,z,delta_x ) 

    % Angular Spetrum Method / Convolution method / 2 fft method
    % without padarray
    % z diffraction distance
    % I0 diffracted image
    % The sampling period and wavelength are in the function

% image input

[N,M]=size(I0); % total sampling number

% parameters setup

% delta_x=10*lambda; %sampling period 10wavelength
% delta_x = 3.8e-3;

delta_y=delta_x;
delta_fx=1/(M*delta_x);% frequency resolution related to sampling periods
delta_fy=1/(N*delta_y);
% z=20; %diffraction distance (mm)

% image processed for diffraction
ift_I0=fftshift(ifft2(fftshift(I0)));

% Spatial Frequency Transfer Function (SFTF)
indice_X=1:M;
indice_Y=1:N;
[indice_x,indice_y]=meshgrid(indice_X,indice_Y);

SFTF=exp(-2i*pi*z.*((1/lambda).^2-((indice_x-M/2-1).*delta_fx).^2-((indice_y-N/2-1).*delta_fy).^2).^0.5);
I_ASM=fftshift(fft2(fftshift(ift_I0.*SFTF)));
end

