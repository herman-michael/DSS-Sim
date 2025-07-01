close all;clear;clc;

load('spectra.mat');

wavelength=spectra.wavelength;
am0=spectra.am0;
sensor_r=spectra.sensor_r;
sensor_g=spectra.sensor_g;
sensor_b=spectra.sensor_b;

am0_norm=am0/max(am0);
sensor_bayer=sensor_r+2*sensor_g+sensor_b;
sensor_norm=sensor_bayer/max(sensor_bayer);
full_spectra=am0_norm.*sensor_norm;
full_norm=full_spectra/max(full_spectra);

wave=[450 full_norm(find(wavelength==450))];
wave=[wave;[wavelength(find(full_norm==1)) full_norm(find(full_norm==1))]];
wave=[wave;[550 full_norm(find(wavelength==550))]];
wave=[wave;[600 full_norm(find(wavelength==600))]];
wave=[wave;[650 full_norm(find(wavelength==650))]];

wave

hold all;
plot(wavelength,am0_norm,'k--');
plot(wavelength,sensor_r,'r');
plot(wavelength,sensor_g,'g');
plot(wavelength,sensor_b,'b');
plot(wavelength,full_norm,'k','Linewidth',2);
plot(wave(:,1),wave(:,2),'ro','Linewidth',2);
xlabel('Wavelength (\lambda)');
ylabel('Normalized Weight');
legend('AM0','R','G','B','Full Response','Wave Samples',...
    'location','SouthEast');