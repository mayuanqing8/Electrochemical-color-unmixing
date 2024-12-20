% open the leica lif file and seperate the channels
clear
file=bfopen();
im_width=size(file{1,1}{1,1},1);
im_high=size(file{1,1}{1,1},2);
total_frame_length=size(file{1,1},1);
green=zeros(im_high,im_width,size(file{1,1},1)/3);
red=zeros(im_high,im_width,size(file{1,1},1)/3);
far_red=zeros(im_high,im_width,size(file{1,1},1)/3);
n=1;
for i=1:3:total_frame_length
far_red(:,:,n)=file{1,1}{i,1};
red(:,:,n)=file{1,1}{i+1,1};
green(:,:,n)=file{1,1}{i+2,1};
n=n+1;
end
plot(squeeze(mean(far_red,[1 2])),'-*');
%% crop a single electrochemcial response cycle in the time lapse.
crop_start=2;
crop_stop=20;
green=green(:,:,crop_start:crop_stop);
red=red(:,:,crop_start:crop_stop);
far_red=far_red(:,:,crop_start:crop_stop);
frame_length=size(green,3);

% drift correction
ref_im=green(:,:,1);
for i=2:size(green,3)
target_im=green(:,:,i);
fft_one=fft2(ref_im);
fft_two=fft2(target_im);
[shift, ~] = dftregistration(fft_one,fft_two,20);
green(:,:,i)=imtranslate(target_im,[shift(4),shift(3)]);
end

ref_im=red(:,:,1);
for i=2:size(red,3)
target_im=red(:,:,i);
fft_one=fft2(ref_im);
fft_two=fft2(target_im);
[shift, ~] = dftregistration(fft_one,fft_two,20);
red(:,:,i)=imtranslate(target_im,[shift(4),shift(3)]);
end

ref_im=far_red(:,:,1);
for i=2:size(far_red,3)
target_im=far_red(:,:,i);
fft_one=fft2(ref_im);
fft_two=fft2(target_im);
[shift, ~] = dftregistration(fft_one,fft_two,20);
far_red(:,:,i)=imtranslate(target_im,[shift(4),shift(3)]);
end

% collect the reference spectral from the composite imgage of SUM and STD
% of the movie.
green_gaus=imgaussfilt3(green,0.1);
green_SUM=mean(green_gaus,3);
green_STD=var(green_gaus,0,3);
green_STD_lo_lim=prctile(green_STD,1,'all');                                                                                                                                                                                             
green_STD_up_lim=prctile(green_STD,99.9,'all');
green_STD(green_STD>green_STD_up_lim)=green_STD_up_lim;
green_STD(green_STD<green_STD_lo_lim)=NaN;
green_SUM_lo_lim=prctile(green_SUM,1,'all');
green_SUM_up_lim=prctile(green_SUM,99.9,'all');
green_SUM(green_SUM>green_SUM_up_lim)=green_SUM_up_lim;
green_SUM(green_SUM<green_SUM_lo_lim)=NaN;
green_diff=diff(green_gaus(:,:,2:14),1,3);
green_diff=sum(green_diff,3);
green_diff_lo_lim=prctile(green_diff,40,'all');                                                                                                                                                                                             
green_diff_up_lim=prctile(green_diff,99.95,'all');
green_diff(green_diff>green_diff_up_lim)=green_diff_up_lim;
green_diff(green_diff<green_diff_lo_lim)=NaN;

red_gaus=imgaussfilt3(red,0.1);
red_SUM=mean(red_gaus,3);
red_STD=var(red_gaus,0,3);
red_STD_lo_lim=prctile(red_STD,10,'all');                                                                                                                                                                                             
red_STD_up_lim=prctile(red_STD,99.95,'all');
red_STD(red_STD>red_STD_up_lim)=red_STD_up_lim;
red_STD(red_STD<red_STD_lo_lim)=NaN;
red_SUM_lo_lim=prctile(red_SUM,10,'all');
red_SUM_up_lim=prctile(red_SUM,99.999,'all');
red_SUM(red_SUM>red_SUM_up_lim)=red_SUM_up_lim;
red_SUM(red_SUM<red_SUM_lo_lim)=NaN;
red_diff=diff(red_gaus(:,:,2:14),1,3);
red_diff=sum(red_diff,3);
red_diff_lo_lim=prctile(red_diff,50,'all');                                                                                                                                                                                             
red_diff_up_lim=prctile(red_diff,99.9999,'all');
red_diff(red_diff>red_diff_up_lim)=red_diff_up_lim;
red_diff(red_diff<red_diff_lo_lim)=NaN;

far_red_gaus=imgaussfilt3(far_red,0.1);
far_red_SUM=mean(far_red_gaus,3);
far_red_STD=var(far_red_gaus,0,3);
far_red_STD_lo_lim=prctile(far_red_STD,10,'all');                                                                                                                                                                                             
far_red_STD_up_lim=prctile(far_red_STD,99.9,'all');
far_red_STD(far_red_STD>far_red_STD_up_lim)=far_red_STD_up_lim;
far_red_STD(far_red_STD<far_red_STD_lo_lim)=NaN;
far_red_SUM_lo_lim=prctile(far_red_SUM,10,'all');
far_red_SUM_up_lim=prctile(far_red_SUM,99.999,'all');
far_red_SUM(far_red_SUM>far_red_SUM_up_lim)=far_red_SUM_up_lim;
far_red_SUM(far_red_SUM<far_red_SUM_lo_lim)=NaN;
far_red_diff=diff(far_red_gaus(:,:,2:14),1,3);
far_red_diff=sum(far_red_diff,3);
far_red_diff_lo_lim=prctile(far_red_diff,40,'all');                                                                                                                                                                                             
far_red_diff_up_lim=prctile(far_red_diff,99.99,'all');
far_red_diff(far_red_diff>far_red_diff_up_lim)=far_red_diff_up_lim;
far_red_diff(far_red_diff<far_red_diff_lo_lim)=NaN;

green_SUM=rescale(green_SUM,0,255);
green_STD=rescale(green_STD,0,255);
green_diff=rescale(green_diff,0,255);
RGB=cat(3, green_SUM,green_STD,green_diff );
figure (1)
imshow(uint8(RGB));
set(gcf, 'Position', get(0, 'Screensize'));
uiwait(msgbox(sprintf('4 click for GFP then 4 for Vimentin')));
[X_cors_G,Y_cors_G] = ginput_white_cursor(8);
close all

red_SUM=rescale(red_SUM,0,255);
red_STD=rescale(red_STD,0,255);
red_diff=rescale(red_diff,0,255);
RGB=cat(3, red_SUM,red_STD,red_diff );
figure (1)
imshow(uint8(RGB));
set(gcf, 'Position', get(0, 'Screensize'));
uiwait(msgbox(sprintf('4 click for mcherry then 4 for Transferin')));
[X_cors_R,Y_cors_R] = ginput_white_cursor(8);
close all

far_red_SUM=rescale(far_red_SUM,0,255);
far_red_STD=rescale(far_red_STD,0,255);
far_red_diff=rescale(far_red_diff,0,255);
RGB=cat(3, far_red_SUM,far_red_STD,far_red_diff );
figure (1)
imshow(uint8(RGB));
set(gcf, 'Position', get(0, 'Screensize'));
uiwait(msgbox(sprintf('4 click for actin then 4 for Pax')));
[X_cors_FR,Y_cors_FR] = ginput_white_cursor(8);
close all

%% unmixing
[GFP,Vimentin,mCherry,Transferin,PAX,actin]=six_color_fun(X_cors_G,X_cors_R,X_cors_FR,Y_cors_G,Y_cors_R,Y_cors_FR,green,red,far_red);

%% save the unmixed image and reference pixel coordinates
save Reference_cor.mat X_cors_G X_cors_R X_cors_FR Y_cors_G Y_cors_R Y_cors_FR;
t = Tiff('GFP.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(GFP));
t.close();

t = Tiff('Vimentin.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(Vimentin));
t.close();

t = Tiff('mCherry.tif','w');
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

t = Tiff('Actin4.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(actin));
t.close();

t = Tiff('PAX4.tif','w');
tagstruct.ImageLength     = size(red,1);
tagstruct.ImageWidth      = size(red,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(PAX));
t.close();