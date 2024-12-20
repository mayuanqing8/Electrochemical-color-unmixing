%% open the experimental file
clear
file=bfopen();
z_plane_num=size(file,1);
im_width=size(file{1,1}{1,1},1);
im_high=size(file{1,1}{1,1},2);
%% cholse the z plane (this is to collect the experimental EC spectra as reference)
z=2;
total_frame_length=size(file{z,1},1);
green=zeros(im_high,im_width,size(file{z,1},1));

for i=1:total_frame_length
green(:,:,i)=file{z,1}{i,1};
end
plot(squeeze(mean(green,[1 2])),'-*');

% crop a single EC cycle in time lapse.
frame_length=size(green,3);
crop_start=12;
crop_stop=50;
green=green(:,:,crop_start:crop_stop);
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

% set up a RGB imagee using the SUM and and STD of the movie for refernce collection. 
green_gaus=imgaussfilt3(green,0.5);
green_SUM=mean(green_gaus,3);
green_STD=var(green_gaus,0,3);
green_STD_lo_lim=prctile(green_STD,1,'all');                                                                                                                                                                                             
green_STD_up_lim=prctile(green_STD,99.99,'all');
green_STD(green_STD>green_STD_up_lim)=green_STD_up_lim;
green_STD(green_STD<green_STD_lo_lim)=NaN;
green_SUM_lo_lim=prctile(green_SUM,1,'all');
green_SUM_up_lim=prctile(green_SUM,99.999,'all');
green_SUM(green_SUM>green_SUM_up_lim)=green_SUM_up_lim;
green_SUM(green_SUM<green_SUM_lo_lim)=NaN;

figure (1)
imagesc(green_SUM);
axis square 
colormap gray
title('raw green sum');

figure (2)
imagesc(green_STD);
colormap gray
axis square
title('raw green STD');
%
close all
green_SUM=rescale(green_SUM,0,255);
green_STD=rescale(green_STD,0,255);
green_diff=zeros(im_high,im_width);
RGB=cat(3, green_STD,green_SUM,green_diff );
figure (1)
imshow(uint8(RGB));
set(gcf, 'Position', get(0, 'Screensize'));
uiwait(msgbox(sprintf('8 click for GFP then 8 for Vimentin')));
[X_cors_G,Y_cors_G] = ginput_white_cursor(16);
close all
%% geting the EC spectra
frame_length=size(green,3);
X_cors_G=round(X_cors_G);
Y_cors_G=round(Y_cors_G);
GFP_ref_1=mean(green(Y_cors_G(1)-3:Y_cors_G(1)+3,X_cors_G(1)-3:X_cors_G(1)+3,1:frame_length),[1 2]);
GFP_ref_2=mean(green(Y_cors_G(2)-3:Y_cors_G(2)+3,X_cors_G(2)-3:X_cors_G(2)+3,1:frame_length),[1 2]);
GFP_ref_3=mean(green(Y_cors_G(3)-3:Y_cors_G(3)+3,X_cors_G(3)-3:X_cors_G(3)+3,1:frame_length),[1 2]);
GFP_ref_4=mean(green(Y_cors_G(4)-3:Y_cors_G(4)+3,X_cors_G(4)-3:X_cors_G(4)+3,1:frame_length),[1 2]);
GFP_ref_5=mean(green(Y_cors_G(5)-3:Y_cors_G(5)+3,X_cors_G(5)-3:X_cors_G(5)+3,1:frame_length),[1 2]);
GFP_ref_6=mean(green(Y_cors_G(6)-3:Y_cors_G(6)+3,X_cors_G(6)-3:X_cors_G(6)+3,1:frame_length),[1 2]);
GFP_ref_7=mean(green(Y_cors_G(7)-3:Y_cors_G(7)+3,X_cors_G(7)-3:X_cors_G(7)+3,1:frame_length),[1 2]);
GFP_ref_8=mean(green(Y_cors_G(8)-3:Y_cors_G(8)+3,X_cors_G(8)-3:X_cors_G(8)+3,1:frame_length),[1 2]);
GFP_ref_com=squeeze(cat(1,GFP_ref_1,GFP_ref_2,GFP_ref_3,GFP_ref_4,GFP_ref_5,GFP_ref_6,GFP_ref_7,GFP_ref_8));

Vimentin_ref_1=mean(green(Y_cors_G(9)-3:Y_cors_G(9)+3,X_cors_G(9)-3:X_cors_G(9)+3,1:frame_length),[1 2]);
Vimentin_ref_2=mean(green(Y_cors_G(10)-3:Y_cors_G(10)+3,X_cors_G(10)-3:X_cors_G(10)+3,1:frame_length),[1 2]);
Vimentin_ref_3=mean(green(Y_cors_G(11)-3:Y_cors_G(11)+3,X_cors_G(11)-3:X_cors_G(11)+3,1:frame_length),[1 2]);
Vimentin_ref_4=mean(green(Y_cors_G(12)-3:Y_cors_G(12)+3,X_cors_G(12)-3:X_cors_G(12)+3,1:frame_length),[1 2]);
Vimentin_ref_5=mean(green(Y_cors_G(13)-3:Y_cors_G(13)+3,X_cors_G(13)-3:X_cors_G(13)+3,1:frame_length),[1 2]);
Vimentin_ref_6=mean(green(Y_cors_G(14)-3:Y_cors_G(14)+3,X_cors_G(14)-3:X_cors_G(14)+3,1:frame_length),[1 2]);
Vimentin_ref_7=mean(green(Y_cors_G(15)-3:Y_cors_G(15)+3,X_cors_G(15)-3:X_cors_G(15)+3,1:frame_length),[1 2]);
Vimentin_ref_8=mean(green(Y_cors_G(16)-3:Y_cors_G(16)+3,X_cors_G(16)-3:X_cors_G(16)+3,1:frame_length),[1 2]);
Vimentin_ref_com=squeeze(cat(1,Vimentin_ref_1,Vimentin_ref_2,Vimentin_ref_3,Vimentin_ref_4,Vimentin_ref_5,...
    Vimentin_ref_6,Vimentin_ref_7,Vimentin_ref_8));

GFP_mean=mean(GFP_ref_com,1);
GFP_smooth=smoothdata(GFP_mean,'movmean',3);
Vimentin_mean=mean(Vimentin_ref_com,1);
Vimentin_smooth=smoothdata(Vimentin_mean,'movmean',3);

figure (4)
plot(GFP_smooth);
hold on
plot(Vimentin_smooth);
hold off
legend('GFP','Vimentin');

%% single pixel simulation using experimental data (varying the ratio of two dyes)
frame_length=numel(Vimentin_smooth);
noise_example=(randn(frame_length,1)*300)'; % add norm distributed noise.
im_size=100;
ref_1=GFP_smooth;
ref_2=Vimentin_smooth;
temporal_signal_combined=cat(2,ref_1',ref_2');
mixture=0.5*ref_1+0.5*ref_2+noise_example; 
% this is simulated signal with nose using realistic experimental data;
options = optimset('TolX',1e-10);
[fitted,ssm,resi]=lsqnonneg(temporal_signal_combined,mixture',options);
% fittted contains the fiting resulted fractions of ref 1 and ref 2.
subplot(2,1,1)
plot(ref_1,'b-','LineWidth',1);
hold on 
plot(ref_2,'r-','LineWidth',1);
plot(mixture,'k*','LineWidth',1,MarkerSize=6);
% plot(resi,'m-.');
plot(fitted(1)*ref_1+fitted(2)*ref_2,'g','LineWidth',1);
hold off
ylabel('Intensity');
legend('Ref-1','Ref-2','Mix','Fitted','Fontsize',11);
legend boxoff
set(gca,'fontsize',11)
subplot(2,1,2)
plot(resi,'-k*','LineWidth',0.75,MarkerSize=6);
ylim([-1000 1000]);
yticks([-1000,0,1000]);
ylabel('Residuals');
xlabel('Frames');
pbaspect([5, 1, 1]);
set(gca,'fontsize',11)
%% spatially varying the fraction of ref 1 and ref 2 in horizontal and 
% vertial direction as sine waves and use fitting to recover the exact fraction at input. 
ground_truth_pattern_1=zeros(im_size,im_size,1);
ground_truth_pattern_2=zeros(im_size,im_size,1);
noise_im=zeros(im_size,im_size,1);
mix_image=zeros(im_size,im_size,frame_length);

for i=1:im_size
    for j=1:im_size
        ground_truth_pattern_1(i,j,:)=(sin(0.3*(j-1))+1)*0.5;
        ground_truth_pattern_2(i,j,:)=(cos(0.3*(1-i)+2)+1)*0.5;
        noise_im(i,j,:)=randn(1)*i*25;
        mix_image(i,j,:)=ref_1*((sin(0.3*(j-1))+1)*0.5)+ref_2*((cos(0.3*(1-i)+2)+1)*0.5)...
            +(randn(frame_length,1)*i*30)';
    end
end

figure (1)
imagesc(ground_truth_pattern_1);
colormap gray
axis square

figure (2)
imagesc(ground_truth_pattern_2);
colormap gray
axis square

figure (3)
imagesc(mean(mix_image,3)); % mix of ground truth ref 1 and 2;
colormap gray
axis square

figure (4)
imagesc(noise_im);
colormap gray
axis square

Fractions_lsq=zeros(im_size,im_size,2);
tic
for i=1:im_size    
    for j=1:im_size
            pixel_temp=squeeze(mix_image(i,j,:));
            [Fractions_lsq(i,j,:),~,~]=lsqnonneg(temporal_signal_combined,double(pixel_temp));
    end
end
toc

ref_1_fitted=rescale(Fractions_lsq(:,:,1),0,255);
ref_2_fitted=rescale(Fractions_lsq(:,:,2),0,255);

% unmixed result
figure (5)
imshowpair(ref_1_fitted,ref_2_fitted,'montage');
title('switching                             non switching');
set(gcf, 'Position', get(0, 'Screensize'));

%% save the images;
ground_truth_pattern_1=rescale(ground_truth_pattern_1,0,255);
ground_truth_pattern_2=rescale(ground_truth_pattern_2,0,255);
noise_im=rescale(noise_im,0,255);

save ref_profile_aquired.mat GFP_ref_com Vimentin_ref_com X_cors_G Y_cors_G
t = Tiff('Fitted_ref_1.tif','w');
tagstruct.ImageLength     = size(mix_image,1);
tagstruct.ImageWidth      = size(mix_image,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(ref_1_fitted));
t.close();

t = Tiff('Fitted_ref_2.tif','w');
tagstruct.ImageLength     = size(mix_image,1);
tagstruct.ImageWidth      = size(mix_image,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(ref_2_fitted));
t.close();


t = Tiff('ground_truth_patteren_1.tif','w');
tagstruct.ImageLength     = size(mix_image,1);
tagstruct.ImageWidth      = size(mix_image,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(ground_truth_pattern_1));
t.close();

t = Tiff('ground_truth_patteren_2.tif','w');
tagstruct.ImageLength     = size(mix_image,1);
tagstruct.ImageWidth      = size(mix_image,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(ground_truth_pattern_2));
t.close();

t = Tiff('noise.tif','w');
tagstruct.ImageLength     = size(mix_image,1);
tagstruct.ImageWidth      = size(mix_image,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct);
t.write(uint16(noise_im));
t.close();
