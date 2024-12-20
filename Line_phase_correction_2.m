function [raw_2]=Line_phase_correction_2(raw_1)
cor_length=4;
im_high=size(raw_1,1);
im_width=size(raw_1,2);
raw_2=zeros(im_high,im_width);
range=-cor_length:cor_length;
for i=1:2:im_high-1
    if sum(raw_1(i,:))>300
        raw_1_odd=raw_1(i,:);
        raw_1_even=raw_1(i+1,:);
        acf=xcorr(raw_1_odd(:),raw_1_even(:),cor_length);
        f=fit(range',double(acf),'gauss1');
        raw_1_even=shift(raw_1_even,[0,f.b1]);
        raw_2(i,:)=raw_1_odd;
        raw_2(i+1,:)=raw_1_even;
    end
end

% plot(f,range,acf);