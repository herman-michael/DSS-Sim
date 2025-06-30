function noisy_image = addNoise(im)
%% Define FullWell and sensor Resolution
fw = 8180; % e-
RNStandardDeviation_e = 3; % e-
bit_level = 8;
%% Load input image and allocate output image
im = double(im);
[h,w] = size(im);
noisy_image = im;
for u = 1:h
    for v = 1:w
            meanSignal_e=graylevel2photoelectrons(im(u,v),fw,bit_level);
            sensorOut_e=custome_poissrnd(meanSignal_e);
            readNoise_e=RNStandardDeviation_e.*randn;
            noisy_image(u,v) = photoelectrons2graylevel(sensorOut_e+...
                readNoise_e,fw,bit_level);
    end
end
noisy_image = uint8(noisy_image);