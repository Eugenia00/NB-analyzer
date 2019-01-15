function out_compare = compare_cellmovement(DataMat_c,rotations)
%%
Mats_R = cell(1);
out_compare = cell(1);
Mats_R(1) = {DataMat_c};
C = 3;
R = double(idivide(size(rotations,2),int8(C),'ceil'));

figure('Position', [300, 400, 1400, 700])

v = var(DataMat_c,0,3);
m = mean(DataMat_c,3);
B = v./m;

e = B-1;
e_plot = e;

maxe = 1;
e_plot(e_plot>(maxe)) = (maxe);
e_plot(e_plot<0) = 0;

a(1) = subplot(R,C,1);
hold on
imshow(e_plot,[],'Colormap',autumn)
colorbar
set(a(1),'Position',get(a(1),'Position'),'Tag','1')

for i = 2:size(rotations,2)
    
    rotation = rotations{i};
    DataMatR = CorrectCellMovement2(DataMat_c,rotation(1),rotation(2),rotation(3),rotation(4));
    Mats_R(i) = {DataMatR};
    
    
    v = var(DataMatR,0,3);
    m = mean(DataMatR,3);
    B = v./m;
    
    e = B-1;
    e_plot = e;
    
    maxe = 1;
    e_plot(e_plot>(maxe)) = (maxe);
    e_plot(e_plot<0) = 0;
    
    
    a(i) = subplot(R,C,i);
    hold on
    imshow(e_plot,[],'Colormap',autumn)
    colorbar
    set(a(i),'Position',get(a(i),'Position'),'Tag',num2str(i))
    
    
end

%beep; pause(0.5); beep; pause(0.5); beep;
load handel.mat;
sound(y(2000:17000), 2*Fs);

waitforbuttonpress;
out_tag = str2double(get(gca,'tag'));
out_compare(1) = Mats_R(out_tag);
out_compare(2) = rotations(out_tag);



%r = idivide(i,int8(C),'ceil');
%c = rem(i,C);


% load handel.mat;
% sound(y(2000:17000), 2*Fs);