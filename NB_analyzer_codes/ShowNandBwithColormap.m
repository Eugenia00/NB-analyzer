function [B,N,e,n] = ShowNandBwithColormap(DataMat)

% % load .lif experiment file
% 
% exp_name = 'C:\Users\cammarota.eugenia\Documents\NumbersAndBrightness\experiments\20170118_3617_NB.lif';
% current_dir = cd;
% cd('C:\Users\cammarota.eugenia\Documents\MATLAB\bfmatlab');
% data = bfopen(exp_name);
% cd(current_dir)
% %%
% 
% % 1 - 2p0AOTF_zoom4_PostDex_0
% % 2 - 0p5AOTF_zoom4_PostDex_1
% % 3 - 0p5AOTF_zoom4_PostDex_2
% % 4 - 0p5AOTF_zoom4_PostDex_3
% % 5 - 0p5AOTF_zoom4_PostDex_4
% % 6 - 0p5AOTF_zoom4_PostDex_5
% % 7 - 0p5AOTF_zoom4_PostDex_6
% % 8 - 0p5AOTF_zoom4_PostDex_7
% % 9 - 0p5AOTF_zoom4_PostDex_8
% % 10 - 0p5AOTF_zoom4_PostDex_9
% % 11 - 0p5AOTF_zoom4_PostDex_10
% % 12 - black_sample
% % 13 - black_nosample
% 
% N_series = 11;
% 
% adress =  data{N_series,1}{1,2};
% ind_semicol = strfind(data{N_series,1}{1,2},';');
% name = adress(ind_semicol(1)+2:ind_semicol(2)-1);
% 
% fprintf('--- Series # %d ---\n',N_series)
% 
% start_from = 7; % excude first 6 frames for fast bleaching
% 
% % load time series in matrix
% T_tot = size(data{N_series,1},1)/2;
% Im = data{N_series,1}{start_from,1};
% DataMat = nan([size(Im),T_tot-start_from+1]);
% DataMat_bf = nan([size(Im),T_tot-start_from+1]);
% 
% for i = start_from:T_tot
%     Im = data{N_series,1}{i*2-1,1}; %fluorescence
%     Im_bf = data{N_series,1}{i*2,1}; %fluorescence
%     DataMat(:,:,i-start_from+1) = Im;
%     DataMat_bf(:,:,i-start_from+1) = Im_bf;
% end
% %
% %DataMat = DataMat2R10;


% intTot = squeeze(mean(mean(DataMat,1),2));
% bleach = (intTot(1)-intTot(end))/intTot(1)*100;
% fprintf('Intensity first to last frame: -%.1f%% \n',bleach)


v = var(DataMat,0,3);
m = mean(DataMat,3);
B = v./m;
N = m.^2./v;

e = B-1;
%n = (m.^2)./(v - m);
n = N.*B./e;


B_plot = B;
e_plot = e;
N_plot = N;
n_plot = n;

maxB = 5;
maxe = 1;
B_plot(B_plot>maxB) = maxB;
e_plot(e_plot>(maxe)) = (maxe);
e_plot(e_plot<0) = 0;
n_plot(n_plot<0) = 0;
n_plot(n_plot>100) = 100;


figure('Position', [300, 400, 1400, 700])

a(1) = subplot(2,3,1); 
%imshow(Im_ave_bf,[])
a(4) = subplot(2,3,4); 
imshow(m,[])

a(2) = subplot(2,3,2);
imshow(B_plot,[0 2],'Colormap',autumn)
colorbar
set(a(2),'Position',get(a(2),'Position'))

a(5) = subplot(2,3,5);
imshow(N_plot,[],'Colormap',autumn)
colorbar
set(a(5),'Position',get(a(5),'Position'))

a(3) = subplot(2,3,3);
imshow(e_plot,[],'Colormap',autumn)
colorbar
set(a(3),'Position',get(a(3),'Position'))

a(6) = subplot(2,3,6);
imshow(n_plot,[],'Colormap',autumn)
colorbar
set(a(6),'Position',get(a(6),'Position'))




shg

end


%% before 27 04 2017
% maxB = 1.2;
% B(B>maxB) = maxB;
% e(e>(maxB-1)) = (maxB-1);
% e(e<0) = 0;
% n(n<0) = 0;
% n(n>100) = 100;
% 
% 
% figure('Position', [300, 400, 1400, 700])
% 
% a(1) = subplot(2,3,1); 
% imshow(Im_ave_bf,[])
% a(4) = subplot(2,3,4); 
% imshow(m,[])
% 
% a(2) = subplot(2,3,2);
% imshow(B,[],'Colormap',autumn)
% colorbar
% set(a(2),'Position',get(a(2),'Position'))
% 
% a(5) = subplot(2,3,5);
% imshow(N,[],'Colormap',autumn)
% colorbar
% set(a(5),'Position',get(a(5),'Position'))
% 
% a(3) = subplot(2,3,3);
% imshow(e,[],'Colormap',autumn)
% colorbar
% set(a(3),'Position',get(a(3),'Position'))
% 
% a(6) = subplot(2,3,6);
% imshow(n,[],'Colormap',autumn)
% colorbar
% set(a(6),'Position',get(a(6),'Position'))