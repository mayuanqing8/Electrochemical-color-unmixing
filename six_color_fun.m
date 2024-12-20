function [GFP,Vimentin,mCherry,Transferin,PAX,actin]=six_color_fun(X_cors_G,X_cors_R,X_cors_FR,Y_cors_G,Y_cors_R,Y_cors_FR,green,red,far_red)
im_width=size(green,1);
im_high=size(green,2);
frame_length=size(green,3);
X_cors_G=round(X_cors_G);
Y_cors_G=round(Y_cors_G);
GFP_ref_1=median(green(Y_cors_G(1)-2:Y_cors_G(1)+2,X_cors_G(1)-2:X_cors_G(1)+2,1:frame_length),[1 2]);
GFP_ref_2=median(green(Y_cors_G(2)-2:Y_cors_G(2)+2,X_cors_G(2)-2:X_cors_G(2)+2,1:frame_length),[1 2]);
GFP_ref_3=median(green(Y_cors_G(3)-2:Y_cors_G(3)+2,X_cors_G(3)-2:X_cors_G(3)+2,1:frame_length),[1 2]);
GFP_ref_4=median(green(Y_cors_G(4)-2:Y_cors_G(4)+2,X_cors_G(4)-2:X_cors_G(4)+2,1:frame_length),[1 2]);
GFP_ref_com=median(cat(1,GFP_ref_1,GFP_ref_2,GFP_ref_3,GFP_ref_4));

Vimentin_ref_1=median(green(Y_cors_G(5)-2:Y_cors_G(5)+2,X_cors_G(5)-2:X_cors_G(5)+2,1:frame_length),[1 2]);
Vimentin_ref_2=median(green(Y_cors_G(6)-2:Y_cors_G(6)+2,X_cors_G(6)-2:X_cors_G(6)+2,1:frame_length),[1 2]);
Vimentin_ref_3=median(green(Y_cors_G(7)-2:Y_cors_G(7)+2,X_cors_G(7)-2:X_cors_G(7)+2,1:frame_length),[1 2]);
Vimentin_ref_4=median(green(Y_cors_G(8)-2:Y_cors_G(8)+2,X_cors_G(8)-2:X_cors_G(8)+2,1:frame_length),[1 2]);
Vimentin_ref_com=median(cat(1,Vimentin_ref_1,Vimentin_ref_2,Vimentin_ref_3,Vimentin_ref_4));

X_cors_R=round(X_cors_R);
Y_cors_R=round(Y_cors_R);
mChery_ref_1=median(red(Y_cors_R(1)-2:Y_cors_R(1)+2,X_cors_R(1)-2:X_cors_R(1)+2,1:frame_length),[1 2]);
mChery_ref_2=median(red(Y_cors_R(2)-2:Y_cors_R(2)+2,X_cors_R(2)-2:X_cors_R(2)+2,1:frame_length),[1 2]);
mChery_ref_3=median(red(Y_cors_R(3)-2:Y_cors_R(3)+2,X_cors_R(3)-2:X_cors_R(3)+2,1:frame_length),[1 2]);
mChery_ref_4=median(red(Y_cors_R(4)-2:Y_cors_R(4)+2,X_cors_R(4)-2:X_cors_R(4)+2,1:frame_length),[1 2]);
mChery_ref_com=median(cat(1,mChery_ref_1,mChery_ref_2,mChery_ref_3,mChery_ref_4));

Transferin_ref_1=median(red(Y_cors_R(5)-2:Y_cors_R(5)+2,X_cors_R(5)-2:X_cors_R(5)+2,1:frame_length),[1 2]);
Transferin_ref_2=median(red(Y_cors_R(6)-2:Y_cors_R(6)+2,X_cors_R(6)-2:X_cors_R(6)+2,1:frame_length),[1 2]);
Transferin_ref_3=median(red(Y_cors_R(7)-2:Y_cors_R(7)+2,X_cors_R(7)-2:X_cors_R(7)+2,1:frame_length),[1 2]);
Transferin_ref_4=median(red(Y_cors_R(8)-2:Y_cors_R(8)+2,X_cors_R(8)-2:X_cors_R(8)+2,1:frame_length),[1 2]);
Transferin_ref_com=median(cat(1,Transferin_ref_1,Transferin_ref_2,Transferin_ref_3,Transferin_ref_4));

X_cors_FR=round(X_cors_FR);
Y_cors_FR=round(Y_cors_FR);
Actin_ref_1=median(far_red(Y_cors_FR(1)-2:Y_cors_FR(1)+2,X_cors_FR(1)-2:X_cors_FR(1)+2,1:frame_length),[1 2]);
Actin_ref_2=median(far_red(Y_cors_FR(2)-2:Y_cors_FR(2)+2,X_cors_FR(2)-2:X_cors_FR(2)+2,1:frame_length),[1 2]);
Actin_ref_3=median(far_red(Y_cors_FR(3)-2:Y_cors_FR(3)+2,X_cors_FR(3)-2:X_cors_FR(3)+2,1:frame_length),[1 2]);
Actin_ref_4=median(far_red(Y_cors_FR(4)-2:Y_cors_FR(4)+2,X_cors_FR(4)-2:X_cors_FR(4)+2,1:frame_length),[1 2]);
Actin_ref_com=median(cat(1,Actin_ref_1,Actin_ref_2,Actin_ref_3,Actin_ref_4));

PAX_ref_1=median(far_red(Y_cors_FR(5)-2:Y_cors_FR(5)+2,X_cors_FR(5)-2:X_cors_FR(5)+2,1:frame_length),[1 2]);
PAX_ref_2=median(far_red(Y_cors_FR(6)-2:Y_cors_FR(6)+2,X_cors_FR(6)-2:X_cors_FR(6)+2,1:frame_length),[1 2]);
PAX_ref_3=median(far_red(Y_cors_FR(7)-2:Y_cors_FR(7)+2,X_cors_FR(7)-2:X_cors_FR(7)+2,1:frame_length),[1 2]);
PAX_ref_4=median(far_red(Y_cors_FR(8)-2:Y_cors_FR(8)+2,X_cors_FR(8)-2:X_cors_FR(8)+2,1:frame_length),[1 2]);
PAX_ref_com=median(cat(1,PAX_ref_1,PAX_ref_2,PAX_ref_3,PAX_ref_4));

GFP_ref_com=squeeze(GFP_ref_com);
Vimentin_ref_com=squeeze(Vimentin_ref_com);
mChery_ref_com=squeeze(mChery_ref_com);
Transferin_ref_com=squeeze(Transferin_ref_com);
Actin_ref_com=squeeze(Actin_ref_com);
PAX_ref_com=squeeze(PAX_ref_com);

Fractions_lsq=zeros(im_width,im_high,2);
ref_combined=cat(2,GFP_ref_com,Vimentin_ref_com);
green_gaus=imgaussfilt3(green);
red_gaus=imgaussfilt3(red);
far_red_gaus=imgaussfilt3(far_red);

options = optimset('TolX',1e-10);
parfor i=1:im_width    
    for j=1:im_high
            pixel_temp=squeeze(green_gaus(i,j,:));
            Fractions_lsq(i,j,:)=lsqnonneg(ref_combined,double(pixel_temp),options);
    end
end
GFP_fitted=Fractions_lsq(:,:,1);
Vimentin_fitted=Fractions_lsq(:,:,2);

Fractions_lsq=zeros(im_width,im_high,2);
ref_combined=cat(2,mChery_ref_com,Transferin_ref_com);
parfor i=1:im_width    
    for j=1:im_high
            pixel_temp=squeeze(red_gaus(i,j,:));
            Fractions_lsq(i,j,:)=lsqnonneg(ref_combined,double(pixel_temp),options);
    end
end
mChery_fitted=Fractions_lsq(:,:,1);
Transferin_fitted=Fractions_lsq(:,:,2);

Fractions_lsq=zeros(im_width,im_high,2);
ref_combined=cat(2,Actin_ref_com,PAX_ref_com);
parfor i=1:im_width    
    for j=1:im_high
            pixel_temp=squeeze(far_red_gaus(i,j,:));
            Fractions_lsq(i,j,:)=lsqnonneg(ref_combined,double(pixel_temp),options);
    end
end
Actin_fitted=Fractions_lsq(:,:,1);
PAX_fitted=Fractions_lsq(:,:,2);

figure (4)
plot(GFP_ref_com);
hold on
plot(Vimentin_ref_com);
plot(mChery_ref_com);
plot(Transferin_ref_com);
plot(Actin_ref_com);
plot(PAX_ref_com);
hold off
legend('GFP','Vimentin','mChery','Transferin','Actin','PAX');

sigma=0.05;
GFP_fitted=imgaussfilt(GFP_fitted,sigma*1);
Vimentin_fitted=imgaussfilt(Vimentin_fitted,sigma*1);
mChery_fitted=imgaussfilt(mChery_fitted,sigma*1);
Transferin_fitted=imgaussfilt(Transferin_fitted,sigma*1);
Actin_fitted=imgaussfilt(Actin_fitted,sigma*1);
PAX_fitted=imgaussfilt(PAX_fitted,sigma*1);

GFP_fitted_fitted_lo=prctile(GFP_fitted,0.001,"all");
GFP_fitted_fitted_up=prctile(GFP_fitted,50,"all");
GFP_fitted(GFP_fitted>GFP_fitted_fitted_up)=GFP_fitted_fitted_up;
GFP_fitted(GFP_fitted<GFP_fitted_fitted_lo)=NaN;

Vimentin_fitted_fitted_lo=prctile(Vimentin_fitted,1,"all");
Vimentin_fitted_fitted_up=prctile(Vimentin_fitted,99.99,"all");
Vimentin_fitted(Vimentin_fitted>Vimentin_fitted_fitted_up)=Vimentin_fitted_fitted_up;
Vimentin_fitted(Vimentin_fitted<Vimentin_fitted_fitted_lo)=NaN;
                                                                
mChery_fitted_fitted_lo=prctile(mChery_fitted,1,"all"); 
mChery_fitted_fitted_up=prctile(mChery_fitted,99.5,"all");
mChery_fitted(mChery_fitted>mChery_fitted_fitted_up)=mChery_fitted_fitted_up;
mChery_fitted(mChery_fitted<mChery_fitted_fitted_lo)=NaN;

Transferin_fitted_fitted_lo=prctile(Transferin_fitted,0.1,"all");
Transferin_fitted_fitted_up=prctile(Transferin_fitted,90,"all");
Transferin_fitted(Transferin_fitted>Transferin_fitted_fitted_up)=Transferin_fitted_fitted_up;
Transferin_fitted(Transferin_fitted<Transferin_fitted_fitted_lo)=NaN;

Actin_fitted_fitted_lo=prctile(Actin_fitted,0.00001,"all");
Actin_fitted_fitted_up=prctile(Actin_fitted,70,"all");
Actin_fitted(Actin_fitted>Actin_fitted_fitted_up)=Actin_fitted_fitted_up;
Actin_fitted(Actin_fitted<Actin_fitted_fitted_lo)=NaN;

PAX_fitted_fitted_lo=prctile(PAX_fitted,20,"all");
PAX_fitted_fitted_up=prctile(PAX_fitted,99.999,"all");
PAX_fitted(PAX_fitted>PAX_fitted_fitted_up)=PAX_fitted_fitted_up;
PAX_fitted(PAX_fitted<PAX_fitted_fitted_lo)=NaN;

green_SUM=mean(green,3);
green_STD=var(green,0,3);
green_STD_lo_lim=prctile(green_STD,0.1,'all');                                                                                                                                                                                             
green_STD_up_lim=prctile(green_STD,99.999,'all');
green_STD(green_STD>green_STD_up_lim)=green_STD_up_lim;
green_STD(green_STD<green_STD_lo_lim)=NaN;
green_SUM_lo_lim=prctile(green_SUM,0.1,'all');
green_SUM_up_lim=prctile(green_SUM,99.999,'all');
green_SUM(green_SUM>green_SUM_up_lim)=green_SUM_up_lim;
green_SUM(green_SUM<green_SUM_lo_lim)=NaN;
red_SUM=mean(red,3);
red_STD=var(red,0,3);
red_STD_lo_lim=prctile(red_STD,0.1,'all');                                                                                                                                                                                             
red_STD_up_lim=prctile(red_STD,99.999,'all');
red_STD(red_STD>red_STD_up_lim)=red_STD_up_lim;
red_STD(red_STD<red_STD_lo_lim)=NaN;
red_SUM_lo_lim=prctile(red_SUM,0.1,'all');
red_SUM_up_lim=prctile(red_SUM,99.999,'all');
red_SUM(red_SUM>red_SUM_up_lim)=red_SUM_up_lim;
red_SUM(red_SUM<red_SUM_lo_lim)=NaN;
far_red_SUM=mean(far_red,3);
far_red_STD=var(far_red,0,3);
far_red_STD_lo_lim=prctile(far_red_STD,0.1,'all');                                                                                                                                                                                             
far_red_STD_up_lim=prctile(far_red_STD,99.999,'all');
far_red_STD(far_red_STD>far_red_STD_up_lim)=far_red_STD_up_lim;
far_red_STD(far_red_STD<far_red_STD_lo_lim)=NaN;
far_red_SUM_lo_lim=prctile(far_red_SUM,0.1,'all');
far_red_SUM_up_lim=prctile(far_red_SUM,99.999,'all');
far_red_SUM(far_red_SUM>far_red_SUM_up_lim)=far_red_SUM_up_lim;
far_red_SUM(far_red_SUM<far_red_SUM_lo_lim)=NaN;

GFP=rescale((green_SUM).*GFP_fitted,0,255);
Vimentin=rescale((green_STD).*Vimentin_fitted,0,255);
mCherry=rescale((red_SUM).*mChery_fitted,0,255);
Transferin=rescale(red_STD.*Transferin_fitted,0,255);
actin=rescale(far_red_SUM.*Actin_fitted,0,255);
PAX=rescale(far_red_SUM.*PAX_fitted,0,255);

figure (1)
imshowpair(GFP,Vimentin,'montage');
title('nGFP                                Vimentin  ');
set(gcf, 'Position', get(0, 'Screensize'));

figure (2)
imshowpair(mCherry,Transferin,'montage');
title('mCherry                                  Transferin   ');
set(gcf, 'Position', get(0, 'Screensize'));

figure (3)
imshowpair(actin,PAX,'montage');
title('Actin                              PAX  ');
set(gcf, 'Position', get(0, 'Screensize'));