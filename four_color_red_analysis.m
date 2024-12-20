%load the leica lif file and drift correct
clear
file=bfopen();
im_width=size(file{1,1}{1,1},1);
im_high=size(file{1,1}{1,1},2);
total_frame_length=size(file{1,1},1);
red=zeros(im_high,im_width,size(file{1,1},1)/2);
n=1;
for i=1:2:total_frame_length
red(:,:,n)=file{1,1}{i,1};
n=n+1;
end
ref_im=red(:,:,1);
for i=2:size(red,3)
target_im=red(:,:,i);
fft_one=fft2(ref_im);
fft_two=fft2(target_im);
[shift, ~] = dftregistration(fft_one,fft_two,20);
red(:,:,i)=imtranslate(target_im,[shift(4),shift(3)]);
end

%% collect the reference EC spectral by click 
% on pixels only contains given structure
crop_start=5;
crop_stop=35;
red_crop=red(:,:,crop_start:crop_stop);
frame_length=size(red,3);
%
red_gaus=imgaussfilt3(red_crop,0.7);
red_SUM=mean(red_gaus,3);
red_STD=std(red_gaus,0,3);
red_STD_lo_lim=prctile(red_STD,0.0001,'all');                                                                                                                                                                                             
red_STD_up_lim=prctile(red_STD,99.99,'all');
red_STD(red_STD>red_STD_up_lim)=red_STD_up_lim;
red_STD(red_STD<red_STD_lo_lim)=NaN;
red_SUM_lo_lim=prctile(red_SUM,10,'all');
red_SUM_up_lim=prctile(red_SUM,99.999,'all');
red_SUM(red_SUM>red_SUM_up_lim)=red_SUM_up_lim;
red_SUM(red_SUM<red_SUM_lo_lim)=NaN;
red_diff=diff(red_gaus(:,:,7:11),1,3);
red_diff=-sum(red_diff,3);
red_diff_lo_lim=prctile(red_diff,60,'all');                                                                                                                                                                                             
red_diff_up_lim=prctile(red_diff,99,'all');
red_diff(red_diff>red_diff_up_lim)=red_diff_up_lim;
red_diff(red_diff<red_diff_lo_lim)=NaN;


close all
red_SUM=rescale(red_SUM,0,275);
red_STD=rescale(red_STD,0,255);
red_diff=rescale(red_diff,0,125);
RGB=cat(3, red_STD, red_SUM,red_diff  );
figure (1)
imshow(uint8(RGB));
set(gcf, 'Position', get(0, 'Screensize'));
uiwait(msgbox(sprintf('5 click each for mito,chery,transferin,actin, background')));
[X_cors,Y_cors] = ginput_white_cursor(25);
close all

%% unmixing
[COX_IV,mCherry,Transferin,Actin]=four_color_red_fun(X_cors,Y_cors,red,file);

%% save the image and reference pixel coordinates.
save Reference_cor.mat X_cors Y_cors;
t = Tiff('COX IV.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(COX_IV));
t.close();


t = Tiff('mCherry_2.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(mCherry));
t.close();


t = Tiff('Transferin.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(Transferin));
t.close();


t = Tiff('Actin_2.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(Actin));
t.close();
