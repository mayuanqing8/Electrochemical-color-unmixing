% load the picoquant exported Tiff STED images
clear
close all
STED_Stack = tiffreadVolume('Tubulin PAX SiRDNA Actin U2OS STED Stack.tif');
frame_length=size(STED_Stack,3);
im_width=size(STED_Stack,1);
im_high=size(STED_Stack,2);

% line drift correction
 for j=1:frame_length
 STED_Stack(:,:,j)=Line_phase_correction_2(STED_Stack(:,:,j));
 end

 % frame drift correction
ref_im=STED_Stack(:,:,1);
for i=2:frame_length
target_im=STED_Stack(:,:,i);
fft_one=fft2(ref_im);
fft_two=fft2(target_im);
[frame_shift, ~] = dftregistration(fft_one,fft_two,30);
STED_Stack(:,:,i)=imtranslate(target_im,[frame_shift(4),frame_shift(3)]);
end

movie_SUM_all=sum(STED_Stack,3);
movie_STD=var(STED_Stack,0,3);
movie_STD_lo_lim=prctile(movie_STD,1,'all');                                                                                                                                                                                             
movie_STD_up_lim=prctile(movie_STD,99.999,'all');
movie_STD(movie_STD>movie_STD_up_lim)=movie_STD_up_lim;
movie_STD(movie_STD<movie_STD_lo_lim)=NaN;
movie_SUM_lo_lim=prctile(movie_SUM_all,1,'all');
movie_SUM_up_lim=prctile(movie_SUM_all,99.99,'all');
movie_SUM_all(movie_SUM_all>movie_SUM_up_lim)=movie_SUM_up_lim;
movie_SUM_all(movie_SUM_all<movie_SUM_lo_lim)=NaN;
movie_diff=diff(STED_Stack(:,:,1:end-1),1,3);
movie_diff=-sum(movie_diff,3);
movie_max_min_lo_lim=prctile(movie_diff,50,'all');
movie_max_min_up_lim=prctile(movie_diff,99.9,'all');
movie_diff(movie_diff>movie_max_min_up_lim)=movie_max_min_up_lim;
movie_diff(movie_diff<movie_max_min_lo_lim)=NaN;

figure (1)
imagesc(movie_SUM_all);
% axis square 
colormap gray
title('raw movie sum');

figure (2)
imagesc(movie_STD);
colormap gray
% axis square
title('raw movie STD');

figure (3)
imagesc(movie_diff);
colormap gray;
% axis square 
%% save the fully line frame drift corrected STED Stack (optional)
for j=1:frame_length
 imwrite(uint16(STED_Stack(:,:,j)),'fully_drift_cor_stack.tiff','WriteMode','append');
end
%% collect the reference EC spectra by click 
% on pixels only contain those structures
movie_SUM_all=rescale(movie_SUM_all,0,195);
movie_STD=rescale(movie_STD,0,255);
movie_diff=rescale(movie_diff,0,155);
Blank=zeros(im_width,im_high);
RGB=cat(3, movie_SUM_all,movie_STD,movie_diff );
imshow(uint8(RGB));

set(gcf, 'Position', get(0, 'Screensize'));
uiwait(msgbox(sprintf('click 5 location each dye Tub/Tom20 PAX SiRDNA Actin')));
[X_cors,Y_cors] = ginput_white_cursor(25);
close all
%% use all frames: good for SirDNA

[Tubulin,PAX,SiRDNA,Actin]=four_color_STED_fun(X_cors,Y_cors,STED_Stack);

%%
save Reference_cor.mat X_cors Y_cors;

t = Tiff('tubulin.tif','w');  
tagstruct.ImageLength     = size(STED_Stack,1);
tagstruct.ImageWidth      = size(STED_Stack,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(Tubulin));
t.close();

t = Tiff('PAX_std.tif','w');
tagstruct.ImageLength     = size(STED_Stack,1);
tagstruct.ImageWidth      = size(STED_Stack,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(PAX));
t.close();

t = Tiff('DNA_std.tif','w');
tagstruct.ImageLength     = size(STED_Stack,1);
tagstruct.ImageWidth      = size(STED_Stack,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(SiRDNA));
t.close();

t = Tiff('Actin_nuk.tif','w');
tagstruct.ImageLength     = size(STED_Stack,1);
tagstruct.ImageWidth      = size(STED_Stack,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(Actin));
t.close();

