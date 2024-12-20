function [COX_IV,mCherry,Transferin,Actin]=four_color_red_fun(X_cors,Y_cors,red,file)
im_width=size(red,1);
im_high=size(red,2);
red_2=zeros(im_high,im_width,size(file{1,1},1)/2);
n=1;
total_frame_length=size(file{1,1},1);
for i=1:2:total_frame_length
red_2(:,:,n)=file{1,1}{i+1,1};
n=n+1;
end
ref_im=red_2(:,:,1);
for i=2:size(red_2,3)
target_im=red_2(:,:,i);
fft_one=fft2(ref_im);
fft_two=fft2(target_im);
[shift, ~] = dftregistration(fft_one,fft_two,20);
red_2(:,:,i)=imtranslate(target_im,[shift(4),shift(3)]);
end

frame_length=size(red,3);
X_cors=round(X_cors);
Y_cors=round(Y_cors);
non_switching_ref_1=sum(red(Y_cors(1)-1:Y_cors(1)+1,X_cors(1)-1:X_cors(1)+1,1:frame_length),[1 2]);
non_switching_ref_2=sum(red(Y_cors(2)-1:Y_cors(2)+1,X_cors(2)-1:X_cors(2)+1,1:frame_length),[1 2]);
non_switching_ref_3=sum(red(Y_cors(3)-1:Y_cors(3)+1,X_cors(3)-1:X_cors(3)+1,1:frame_length),[1 2]);
non_switching_ref_4=sum(red(Y_cors(4)-1:Y_cors(4)+1,X_cors(4)-1:X_cors(4)+1,1:frame_length),[1 2]);
non_switching_ref_5=sum(red(Y_cors(5)-1:Y_cors(5)+1,X_cors(5)-1:X_cors(5)+1,1:frame_length),[1 2]);
non_switching_ref_com=mean(cat(1,non_switching_ref_1,non_switching_ref_2,non_switching_ref_3,...
    non_switching_ref_4,non_switching_ref_5));

switching_ref1_1=sum(red(Y_cors(6)-1:Y_cors(6)+1,X_cors(6)-1:X_cors(6)+1,1:frame_length),[1 2]);
switching_ref1_2=sum(red(Y_cors(7)-1:Y_cors(7)+1,X_cors(7)-1:X_cors(7)+1,1:frame_length),[1 2]);
switching_ref1_3=sum(red(Y_cors(8)-1:Y_cors(8)+1,X_cors(8)-1:X_cors(8)+1,1:frame_length),[1 2]);
switching_ref1_4=sum(red(Y_cors(9)-1:Y_cors(9)+1,X_cors(9)-1:X_cors(9)+1,1:frame_length),[1 2]);
switching_ref1_5=sum(red(Y_cors(10)-1:Y_cors(10)+1,X_cors(10)-1:X_cors(10)+1,1:frame_length),[1 2]);

switching_ref1_com=mean(cat(1,switching_ref1_1,switching_ref1_2,switching_ref1_3,switching_ref1_4,...
    switching_ref1_5));

switching_ref2_1=sum(red(Y_cors(11)-1:Y_cors(11)+1,X_cors(11)-1:X_cors(11)+1,1:frame_length),[1 2]);
switching_ref2_2=sum(red(Y_cors(12)-1:Y_cors(12)+1,X_cors(12)-1:X_cors(12)+1,1:frame_length),[1 2]);
switching_ref2_3=sum(red(Y_cors(13)-1:Y_cors(13)+1,X_cors(13)-1:X_cors(13)+1,1:frame_length),[1 2]);
switching_ref2_4=sum(red(Y_cors(14)-1:Y_cors(14)+1,X_cors(14)-1:X_cors(14)+1,1:frame_length),[1 2]);
switching_ref2_5=sum(red(Y_cors(15)-1:Y_cors(15)+1,X_cors(15)-1:X_cors(15)+1,1:frame_length),[1 2]);
switching_ref2_com=mean(cat(1,switching_ref2_1,switching_ref2_2,switching_ref2_3,...
    switching_ref2_4, switching_ref2_5));

switching_ref3_1=sum(red(Y_cors(16)-1:Y_cors(16)+1,X_cors(16)-1:X_cors(16)+1,1:frame_length),[1 2]);
switching_ref3_2=sum(red(Y_cors(17)-1:Y_cors(17)+1,X_cors(17)-1:X_cors(17)+1,1:frame_length),[1 2]);
switching_ref3_3=sum(red(Y_cors(18)-1:Y_cors(18)+1,X_cors(18)-1:X_cors(18)+1,1:frame_length),[1 2]);
switching_ref3_4=sum(red(Y_cors(19)-1:Y_cors(19)+1,X_cors(19)-1:X_cors(19)+1,1:frame_length),[1 2]);
switching_ref3_5=sum(red(Y_cors(20)-1:Y_cors(20)+1,X_cors(20)-1:X_cors(20)+1,1:frame_length),[1 2]);
switching_ref3_com=mean(cat(1,switching_ref3_1,switching_ref3_2,switching_ref3_3,...
    switching_ref3_4,switching_ref3_5));

switching_ref4_1=sum(red(Y_cors(21)-1:Y_cors(21)+1,X_cors(21)-1:X_cors(21)+1,1:frame_length),[1 2]);
switching_ref4_2=sum(red(Y_cors(22)-1:Y_cors(22)+1,X_cors(22)-1:X_cors(22)+1,1:frame_length),[1 2]);
switching_ref4_3=sum(red(Y_cors(23)-1:Y_cors(23)+1,X_cors(23)-1:X_cors(23)+1,1:frame_length),[1 2]);
switching_ref4_4=sum(red(Y_cors(24)-1:Y_cors(24)+1,X_cors(24)-1:X_cors(24)+1,1:frame_length),[1 2]);
switching_ref4_5=sum(red(Y_cors(25)-1:Y_cors(25)+1,X_cors(25)-1:X_cors(25)+1,1:frame_length),[1 2]);
switching_ref4_com=mean(cat(1,switching_ref4_1,switching_ref4_2,switching_ref4_3,...
    switching_ref4_4,switching_ref4_5));

non_switching_ref_com=squeeze(non_switching_ref_com);
switching_ref1_com=squeeze(switching_ref1_com);
switching_ref2_com=squeeze(switching_ref2_com);
switching_ref3_com=squeeze(switching_ref3_com);
switching_ref4_com=squeeze(switching_ref4_com);

Fractions_lsq=zeros(im_width,im_high,5);
non_swtiching_avg_all=non_switching_ref_com./sum(non_switching_ref_com);
swtiching_avg_all_1=switching_ref1_com./sum(switching_ref1_com);
swtiching_avg_all_2=switching_ref2_com./sum(switching_ref2_com);
swtiching_avg_all_3=switching_ref3_com./sum(switching_ref3_com);
swtiching_avg_all_4=switching_ref4_com./sum(switching_ref4_com);

temporal_signal_combined=cat(2,non_swtiching_avg_all,swtiching_avg_all_1,...
    swtiching_avg_all_2,swtiching_avg_all_3,swtiching_avg_all_4);
temporal_signal_combined=double(temporal_signal_combined);

options = optimset('TolX',1e-10);
parfor i=1:im_width    
    for j=1:im_high
            pixel_temp=squeeze(red(j,i,:));
            pixel_temp=pixel_temp./sum(pixel_temp);
            Fractions_lsq(i,j,:)=lsqnonneg(temporal_signal_combined,double(pixel_temp),options);
    end
end

red_STD=var(red,0,3);

switching_fitted_2=Fractions_lsq(:,:,3)';

sigma=0.05;
switching_fitted_2=imgaussfilt(switching_fitted_2,sigma*1);

red_STD_lo_lim=prctile(red_STD,0.0001,'all');                                                                                                                                                                                             
red_STD_up_lim=prctile(red_STD,99.999,'all');
red_STD(red_STD>red_STD_up_lim)=red_STD_up_lim;
red_STD(red_STD<red_STD_lo_lim)=0;

switching_fitted_2(switching_fitted_2==0)=NaN;
                                                                                                                                                                                                                                                                                                                                                                                                                                                      
switching_fitted_2_lo=prctile(switching_fitted_2,0.1,"all");
switching_fitted_2_up=prctile(switching_fitted_2,99.9,"all");
switching_fitted_2(switching_fitted_2>switching_fitted_2_up)=switching_fitted_2_up;
switching_fitted_2(switching_fitted_2<switching_fitted_2_lo)=NaN;

Transferin=rescale((red_STD).*switching_fitted_2,0,255);

red=red+red_2;
frame_length=size(red,3);
X_cors=round(X_cors);
Y_cors=round(Y_cors);
non_switching_ref_1=sum(red(Y_cors(1)-1:Y_cors(1)+1,X_cors(1)-1:X_cors(1)+1,1:frame_length),[1 2]);
non_switching_ref_2=sum(red(Y_cors(2)-1:Y_cors(2)+1,X_cors(2)-1:X_cors(2)+1,1:frame_length),[1 2]);
non_switching_ref_3=sum(red(Y_cors(3)-1:Y_cors(3)+1,X_cors(3)-1:X_cors(3)+1,1:frame_length),[1 2]);
non_switching_ref_4=sum(red(Y_cors(4)-1:Y_cors(4)+1,X_cors(4)-1:X_cors(4)+1,1:frame_length),[1 2]);
non_switching_ref_5=sum(red(Y_cors(5)-1:Y_cors(5)+1,X_cors(5)-1:X_cors(5)+1,1:frame_length),[1 2]);
non_switching_ref_com=mean(cat(1,non_switching_ref_1,non_switching_ref_2,non_switching_ref_3,...
    non_switching_ref_4,non_switching_ref_5));

switching_ref1_1=sum(red(Y_cors(6)-1:Y_cors(6)+1,X_cors(6)-1:X_cors(6)+1,1:frame_length),[1 2]);
switching_ref1_2=sum(red(Y_cors(7)-1:Y_cors(7)+1,X_cors(7)-1:X_cors(7)+1,1:frame_length),[1 2]);
switching_ref1_3=sum(red(Y_cors(8)-1:Y_cors(8)+1,X_cors(8)-1:X_cors(8)+1,1:frame_length),[1 2]);
switching_ref1_4=sum(red(Y_cors(9)-1:Y_cors(9)+1,X_cors(9)-1:X_cors(9)+1,1:frame_length),[1 2]);
switching_ref1_5=sum(red(Y_cors(10)-1:Y_cors(10)+1,X_cors(10)-1:X_cors(10)+1,1:frame_length),[1 2]);

switching_ref1_com=mean(cat(1,switching_ref1_1,switching_ref1_2,switching_ref1_3,switching_ref1_4,...
    switching_ref1_5));

switching_ref2_1=sum(red(Y_cors(11)-1:Y_cors(11)+1,X_cors(11)-1:X_cors(11)+1,1:frame_length),[1 2]);
switching_ref2_2=sum(red(Y_cors(12)-1:Y_cors(12)+1,X_cors(12)-1:X_cors(12)+1,1:frame_length),[1 2]);
switching_ref2_3=sum(red(Y_cors(13)-1:Y_cors(13)+1,X_cors(13)-1:X_cors(13)+1,1:frame_length),[1 2]);
switching_ref2_4=sum(red(Y_cors(14)-1:Y_cors(14)+1,X_cors(14)-1:X_cors(14)+1,1:frame_length),[1 2]);
switching_ref2_5=sum(red(Y_cors(15)-1:Y_cors(15)+1,X_cors(15)-1:X_cors(15)+1,1:frame_length),[1 2]);
switching_ref2_com=mean(cat(1,switching_ref2_1,switching_ref2_2,switching_ref2_3,...
    switching_ref2_4, switching_ref2_5));

switching_ref3_1=sum(red(Y_cors(16)-1:Y_cors(16)+1,X_cors(16)-1:X_cors(16)+1,1:frame_length),[1 2]);
switching_ref3_2=sum(red(Y_cors(17)-1:Y_cors(17)+1,X_cors(17)-1:X_cors(17)+1,1:frame_length),[1 2]);
switching_ref3_3=sum(red(Y_cors(18)-1:Y_cors(18)+1,X_cors(18)-1:X_cors(18)+1,1:frame_length),[1 2]);
switching_ref3_4=sum(red(Y_cors(19)-1:Y_cors(19)+1,X_cors(19)-1:X_cors(19)+1,1:frame_length),[1 2]);
switching_ref3_5=sum(red(Y_cors(20)-1:Y_cors(20)+1,X_cors(20)-1:X_cors(20)+1,1:frame_length),[1 2]);
switching_ref3_com=mean(cat(1,switching_ref3_1,switching_ref3_2,switching_ref3_3,...
    switching_ref3_4,switching_ref3_5));

switching_ref4_1=sum(red(Y_cors(21)-1:Y_cors(21)+1,X_cors(21)-1:X_cors(21)+1,1:frame_length),[1 2]);
switching_ref4_2=sum(red(Y_cors(22)-1:Y_cors(22)+1,X_cors(22)-1:X_cors(22)+1,1:frame_length),[1 2]);
switching_ref4_3=sum(red(Y_cors(23)-1:Y_cors(23)+1,X_cors(23)-1:X_cors(23)+1,1:frame_length),[1 2]);
switching_ref4_4=sum(red(Y_cors(24)-1:Y_cors(24)+1,X_cors(24)-1:X_cors(24)+1,1:frame_length),[1 2]);
switching_ref4_5=sum(red(Y_cors(25)-1:Y_cors(25)+1,X_cors(25)-1:X_cors(25)+1,1:frame_length),[1 2]);
switching_ref4_com=mean(cat(1,switching_ref4_1,switching_ref4_2,switching_ref4_3,...
    switching_ref4_4,switching_ref4_5));

non_switching_ref_com=squeeze(non_switching_ref_com);
switching_ref1_com=squeeze(switching_ref1_com);
switching_ref2_com=squeeze(switching_ref2_com);
switching_ref3_com=squeeze(switching_ref3_com);
switching_ref4_com=squeeze(switching_ref4_com);

Fractions_lsq=zeros(im_width,im_high,5);
non_swtiching_avg_all=non_switching_ref_com./sum(non_switching_ref_com);
swtiching_avg_all_1=switching_ref1_com./sum(switching_ref1_com);
swtiching_avg_all_2=switching_ref2_com./sum(switching_ref2_com);
swtiching_avg_all_3=switching_ref3_com./sum(switching_ref3_com);
swtiching_avg_all_4=switching_ref4_com./sum(switching_ref4_com);

temporal_signal_combined=cat(2,non_swtiching_avg_all,swtiching_avg_all_1,...
    swtiching_avg_all_2,swtiching_avg_all_3,swtiching_avg_all_4);
temporal_signal_combined=double(temporal_signal_combined);

options = optimset('TolX',1e-10);
parfor i=1:im_width    
    for j=1:im_high
            pixel_temp=squeeze(red(j,i,:));
            pixel_temp=pixel_temp./sum(pixel_temp);
            Fractions_lsq(i,j,:)=lsqnonneg(temporal_signal_combined,double(pixel_temp),options);
    end
end

red_SUM=mean(red,3);
red_STD=var(red,0,3);

non_switching_fitted=Fractions_lsq(:,:,1)';
switching_fitted_1=Fractions_lsq(:,:,2)';
switching_fitted_3=Fractions_lsq(:,:,4)';

sigma=0.05;
non_switching_fitted=imgaussfilt(non_switching_fitted,sigma);
switching_fitted_1=imgaussfilt(switching_fitted_1,sigma*10);
switching_fitted_3=imgaussfilt(switching_fitted_3,sigma*1);

red_STD_lo_lim=prctile(red_STD,0.0001,'all');                                                                                                                                                                                             
red_STD_up_lim=prctile(red_STD,99.999,'all');
red_STD(red_STD>red_STD_up_lim)=red_STD_up_lim;
red_STD(red_STD<red_STD_lo_lim)=0;

red_SUM_lo_lim=prctile(red_SUM,0.0001,'all');
red_SUM_up_lim=prctile(red_SUM,99.999,'all');
red_SUM(red_SUM>red_SUM_up_lim)=red_SUM_up_lim;
red_SUM(red_SUM<red_SUM_lo_lim)=0;

non_switching_fitted(non_switching_fitted==0)=NaN;
switching_fitted_1(switching_fitted_1==0)=NaN;
switching_fitted_3(switching_fitted_3==0)=NaN;

non_switching_fitted_lo=prctile(non_switching_fitted,1,"all");
non_switching_fitted_up=prctile(non_switching_fitted,99.9,"all");
non_switching_fitted(non_switching_fitted>non_switching_fitted_up)=non_switching_fitted_up;
non_switching_fitted(non_switching_fitted<non_switching_fitted_lo)=NaN;

switching_fitted_1_lo=prctile(switching_fitted_1,1,"all");
switching_fitted_1_up=prctile(switching_fitted_1,99.99,"all");
switching_fitted_1(switching_fitted_1>switching_fitted_1_up)=switching_fitted_1_up;
switching_fitted_1(switching_fitted_1<switching_fitted_1_lo)=NaN;
                                                                                                                                                                                                                                                                                                                                                                                                                                                            
switching_fitted_3_lo=prctile(switching_fitted_3,0.001,"all");
switching_fitted_3_up=prctile(switching_fitted_3,99.5,"all");
switching_fitted_3(switching_fitted_3>switching_fitted_3_up)=switching_fitted_3_up;
switching_fitted_3(switching_fitted_3<switching_fitted_3_lo)=NaN;

COX_IV=rescale((red_SUM).*non_switching_fitted,0,255);
mCherry=rescale(red_SUM.*switching_fitted_1,0,255);
Actin=rescale(red_SUM.*switching_fitted_3,0,255);


RGB=cat(3, COX_IV,Transferin, mCherry );
Megenta=cat(3,0.7*ones(im_high), 0.7*ones(im_high), 0.1*ones(im_high));

figure (3)
imshowpair(COX_IV,mCherry,'montage');
title('mito-AF594                                 mCherry-H2B   ');
set(gcf, 'Position', get(0, 'Screensize'));

figure (4)
imshowpair(Transferin,Actin,'montage');
title('Transferin-AF568                                 Actin-AF555   ');
set(gcf, 'Position', get(0, 'Screensize'));

overlay_image=figure(2);
imshow(uint8(RGB));
set(overlay_image,'Position',[750,-50,800,800]);
hold on
h = imshow(Megenta);
hold off
set(h, 'AlphaData', Actin*0.005);
axis square;