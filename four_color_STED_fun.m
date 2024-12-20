function [Tubulin,PAX,SiRDNA,Actin]=four_color_STED_fun(X_cors,Y_cors,STED_Stack)
im_width=size(STED_Stack,1);
im_high=size(STED_Stack,2);
frame_length=size(STED_Stack,3)-1;
X_cors=round(X_cors);
Y_cors=round(Y_cors);
non_switching_ref_1=median(STED_Stack(Y_cors(1)-3:Y_cors(1)+3,X_cors(1)-3:X_cors(1)+3,1:frame_length),[1 2]);
non_switching_ref_2=median(STED_Stack(Y_cors(2)-3:Y_cors(2)+3,X_cors(2)-3:X_cors(2)+3,1:frame_length),[1 2]);
non_switching_ref_3=median(STED_Stack(Y_cors(3)-3:Y_cors(3)+3,X_cors(3)-3:X_cors(3)+3,1:frame_length),[1 2]);
non_switching_ref_4=median(STED_Stack(Y_cors(4)-3:Y_cors(4)+3,X_cors(4)-3:X_cors(4)+3,1:frame_length),[1 2]);
non_switching_ref_5=median(STED_Stack(Y_cors(5)-3:Y_cors(5)+3,X_cors(5)-3:X_cors(5)+3,1:frame_length),[1 2]);
non_switching_ref_com=median(cat(1,non_switching_ref_1,non_switching_ref_2,non_switching_ref_3,...
    non_switching_ref_4,non_switching_ref_5));

switching_ref1_1=median(STED_Stack(Y_cors(6)-3:Y_cors(6)+3,X_cors(6)-3:X_cors(6)+3,1:frame_length),[1 2]);
switching_ref1_2=median(STED_Stack(Y_cors(7)-3:Y_cors(7)+3,X_cors(7)-3:X_cors(7)+3,1:frame_length),[1 2]);
switching_ref1_3=median(STED_Stack(Y_cors(8)-3:Y_cors(8)+3,X_cors(8)-3:X_cors(8)+3,1:frame_length),[1 2]);
switching_ref1_4=median(STED_Stack(Y_cors(9)-3:Y_cors(9)+3,X_cors(9)-3:X_cors(9)+3,1:frame_length),[1 2]);
switching_ref1_5=median(STED_Stack(Y_cors(10)-3:Y_cors(10)+3,X_cors(10)-3:X_cors(10)+3,1:frame_length),[1 2]);

switching_ref1_com=median(cat(1,switching_ref1_1,switching_ref1_2,switching_ref1_3,switching_ref1_4,...
    switching_ref1_5));

switching_ref2_1=sum(STED_Stack(Y_cors(11)-3:Y_cors(11)+3,X_cors(11)-3:X_cors(11)+3,1:frame_length),[1 2]);
switching_ref2_2=sum(STED_Stack(Y_cors(12)-3:Y_cors(12)+3,X_cors(12)-3:X_cors(12)+3,1:frame_length),[1 2]);
switching_ref2_3=sum(STED_Stack(Y_cors(13)-3:Y_cors(13)+3,X_cors(13)-3:X_cors(13)+3,1:frame_length),[1 2]);
switching_ref2_4=sum(STED_Stack(Y_cors(14)-3:Y_cors(14)+3,X_cors(14)-3:X_cors(14)+3,1:frame_length),[1 2]);
switching_ref2_5=sum(STED_Stack(Y_cors(15)-3:Y_cors(15)+3,X_cors(15)-3:X_cors(15)+3,1:frame_length),[1 2]);
switching_ref2_com=sum(cat(1,switching_ref2_1,switching_ref2_2,switching_ref2_3,...
    switching_ref2_4, switching_ref2_5));


switching_ref3_1=median(STED_Stack(Y_cors(16)-3:Y_cors(16)+3,X_cors(16)-3:X_cors(16)+3,1:frame_length),[1 2]);
switching_ref3_2=median(STED_Stack(Y_cors(17)-3:Y_cors(17)+3,X_cors(17)-3:X_cors(17)+3,1:frame_length),[1 2]);
switching_ref3_3=median(STED_Stack(Y_cors(18)-3:Y_cors(18)+3,X_cors(18)-3:X_cors(18)+3,1:frame_length),[1 2]);
switching_ref3_4=median(STED_Stack(Y_cors(19)-3:Y_cors(19)+3,X_cors(19)-3:X_cors(19)+3,1:frame_length),[1 2]);
switching_ref3_5=median(STED_Stack(Y_cors(20)-3:Y_cors(20)+3,X_cors(20)-3:X_cors(20)+3,1:frame_length),[1 2]);
switching_ref3_com=median(cat(1,switching_ref3_1,switching_ref3_2,switching_ref3_3,...
    switching_ref3_4,switching_ref3_5));

switching_ref4_1=median(STED_Stack(Y_cors(21)-3:Y_cors(21)+3,X_cors(21)-3:X_cors(21)+3,1:frame_length),[1 2]);
switching_ref4_2=median(STED_Stack(Y_cors(22)-3:Y_cors(22)+3,X_cors(22)-3:X_cors(22)+3,1:frame_length),[1 2]);
switching_ref4_3=median(STED_Stack(Y_cors(23)-3:Y_cors(23)+3,X_cors(23)-3:X_cors(23)+3,1:frame_length),[1 2]);
switching_ref4_4=median(STED_Stack(Y_cors(24)-3:Y_cors(24)+3,X_cors(24)-3:X_cors(24)+3,1:frame_length),[1 2]);
switching_ref4_5=median(STED_Stack(Y_cors(25)-3:Y_cors(25)+3,X_cors(25)-3:X_cors(25)+3,1:frame_length),[1 2]);
switching_ref4_com=median(cat(1,switching_ref4_1,switching_ref4_2,switching_ref4_3,...
    switching_ref4_4,switching_ref4_5));

non_switching_ref_com=squeeze(non_switching_ref_com);
switching_ref1_com=squeeze(switching_ref1_com);
switching_ref2_com=squeeze(switching_ref2_com);
switching_ref3_com=squeeze(switching_ref3_com);
switching_ref4_com=squeeze(switching_ref4_com);

Fractions_lsq=zeros(im_width,im_high,5);
temporal_signal_combined=cat(2,non_switching_ref_com,switching_ref1_com,...
    switching_ref2_com,switching_ref3_com,switching_ref4_com);
temporal_signal_combined=double(temporal_signal_combined);

options = optimset('TolX',1e-10);
parfor i=1:im_width    
    for j=1:im_high
            pixel_temp=squeeze(STED_Stack(i,j,1:frame_length));
            Fractions_lsq(i,j,:)=lsqnonneg(temporal_signal_combined,double(pixel_temp),options);
    end
end

movie_SUM=mean(STED_Stack,3);
movie_STD=std(double(STED_Stack),0,3);

non_switching_fitted=Fractions_lsq(:,:,1);
switching_fitted_1=Fractions_lsq(:,:,2);
switching_fitted_2=Fractions_lsq(:,:,3);
switching_fitted_3=Fractions_lsq(:,:,4);

sigma=0.05;
non_switching_fitted=imgaussfilt(non_switching_fitted,sigma*3);
switching_fitted_1=imgaussfilt(switching_fitted_1,sigma*3);
switching_fitted_2=imgaussfilt(switching_fitted_2,sigma*3);
switching_fitted_3=imgaussfilt(switching_fitted_3,sigma*1);

movie_STD_lo_lim=prctile(movie_STD,0.01,'all');                                                                                                                                                                                             
movie_STD_up_lim=prctile(movie_STD,99.9999,'all');
movie_STD(movie_STD>movie_STD_up_lim)=movie_STD_up_lim;
movie_STD(movie_STD<movie_STD_lo_lim)=0;

movie_SUM_lo_lim=prctile(movie_SUM,0.01,'all');
movie_SUM_up_lim=prctile(movie_SUM,99.9999,'all');
movie_SUM(movie_SUM>movie_SUM_up_lim)=movie_SUM_up_lim;
movie_SUM(movie_SUM<movie_SUM_lo_lim)=0;

non_switching_fitted(non_switching_fitted==0)=NaN;
switching_fitted_1(switching_fitted_1==0)=NaN;
switching_fitted_2(switching_fitted_2==0)=NaN;
switching_fitted_3(switching_fitted_3==0)=NaN;

non_switching_fitted_lo=prctile(non_switching_fitted,0.00001,"all");
non_switching_fitted_up=prctile(non_switching_fitted,80,"all");
non_switching_fitted(non_switching_fitted>non_switching_fitted_up)=non_switching_fitted_up;
non_switching_fitted(non_switching_fitted<non_switching_fitted_lo)=NaN;

switching_fitted_1_lo=prctile(switching_fitted_1,0.00001,"all");
switching_fitted_1_up=prctile(switching_fitted_1,99,"all");
switching_fitted_1(switching_fitted_1>switching_fitted_1_up)=switching_fitted_1_up;
switching_fitted_1(switching_fitted_1<switching_fitted_1_lo)=NaN;
                                                                                                                                                                                                                                                                                                                                                                                                                                                            
switching_fitted_2_lo=prctile(switching_fitted_2,0.001,"all");
switching_fitted_2_up=prctile(switching_fitted_2,99.9,"all");
switching_fitted_2(switching_fitted_2>switching_fitted_2_up)=switching_fitted_2_up;
switching_fitted_2(switching_fitted_2<switching_fitted_2_lo)=NaN;


switching_fitted_3_lo=prctile(switching_fitted_3,0.00001,"all");
switching_fitted_3_up=prctile(switching_fitted_3,1,"all");
switching_fitted_3(switching_fitted_3>switching_fitted_3_up)=switching_fitted_3_up;
switching_fitted_3(switching_fitted_3<switching_fitted_3_lo)=NaN;

Tubulin=rescale(movie_SUM.*non_switching_fitted,0,255);
PAX=rescale(movie_SUM.*switching_fitted_1,0,255);
SiRDNA=rescale(movie_SUM.*switching_fitted_2,0,255);
Actin=rescale(movie_STD.*switching_fitted_3,0,255);
RGB=cat(3, Tubulin,SiRDNA, PAX );
Megenta=cat(3,0.8*ones(im_width,im_high), 0.2*ones(im_width,im_high), 0.8*ones(im_width,im_high));


figure (3)
imshowpair(Tubulin,PAX,'montage');
title('Tubulin-Starred                                  Pax-Atto647N   ');
set(gcf, 'Position', get(0, 'Screensize'));

figure (4)
imshowpair(SiRDNA,Actin,'montage');
title('SiR DNA                                 Actin-Atto655  ');
set(gcf, 'Position', get(0, 'Screensize'));


overlay_image=figure(2);
imshow(uint8(RGB));
set(overlay_image,'Position',[750,-50,800,800]);
hold on
h = imshow(Megenta);
hold off
set(h, 'AlphaData', Actin*0.003);
end
