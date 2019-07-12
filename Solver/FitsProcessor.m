%% Read info from fits files
%% siddharth@deore.in
%dirpath='C:\Users\Siddharth Deore\Desktop\Processed Fits';
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_26_30.5_00020110.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_27_0.5_00020111.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_27_30.5_00020112.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_28_1.5_00020113.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_28_30.5_00020114.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_29_0.5_00020115.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_29_30.5_00020116.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_30_0.5_00020117.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_30_30.5_00020118.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_32_3.5_00020121.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_32_31.5_00020122.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_33_0.5_00020123.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_21_03_33_31.5_00020124.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_02_58_30.5_00020745.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_02_59_0.5_00020746.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_02_59_31.5_00020747.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_00_0.5_00020748.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_00_30.5_00020749.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_00_30.5_00020754.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_01_0.5_00020750.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_01_31.5_00020751.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_02_1.5_00020752.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_02_31.5_00020753.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_03_30.5_00020755.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_04_1.5_00020756.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_22_03_04_31.5_00020757.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_31_31.5_00020947.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_32_0.5_00020948.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_32_30.5_00020949.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_33_0.5_00020950.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_33_31.5_00020951.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_34_1.5_00020952.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_34_30.5_00020953.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_02_35_1.5_00020954.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_10_30.5_00020970.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_11_30.5_00020972.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_12_0.5_00020973.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_12_30.5_00020974.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_13_1.5_00020975.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_13_31.5_00020976.fits');
getData('C:\Users\Siddharth Deore\Desktop\ProcessedFits\27422_2019_03_23_04_14_0.5_00020977.fits');


function [lst,lat,lon,HA,utc,local]= getData(fit)
info=fitsinfo(fit);
c=info.PrimaryData.Keywords;
i1=find(strcmp(c, 'LST'));
i2=find(strcmp(c, 'SITELAT'));
i3=find(strcmp(c, 'SITELONG'));
i4=find(strcmp(c, 'TELEHA'));
i5=find(strcmp(c, 'DATE-OBS'));
i6=find(strcmp(c, 'LOCALTIM'));

lst=cell2mat(c(i1,2));
lat=cell2mat(c(i2,2));
lon=cell2mat(c(i3,2));
HA=cell2mat(c(i4,2));
utc=cell2mat(c(i5,2));
local=cell2mat(c(i6,2));
%fprintf(strcat(lst,', ',lat,', ',lon,', ',HA,', ',utc,', ',local,'\n'));
fprintf(strcat(utc,'\n'));
end