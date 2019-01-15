%% header store data
res = cell(1,1);
res(1,1) = {'lif file'};
res(1,2) = {'type'};
res(1,3) = {'name series'};
res(1,4) = {'N'};
res(1,5) = {'B'};
res(1,6) = {'n'};
res(1,7) = {'e'};
res(1,8) = {'N std'};
res(1,9) = {'B std'};
res(1,10) = {'n std'};
res(1,11) = {'e std'};
res(1,12) = {'bleaching'};
res(1,13) = {'rect cell'};
res(1,14) = {'rect 5 regions'};
res(1,15) = {'movement correction Tmax,Rmax,verifyEvery,everage on'};
res(1,16) = {'Intensity fluo'};
res(1,17) = {'Intensity bf'};
res(1,18) = {'n Boxcar'};
res(1,19) = {'e Boxcar'};

%% load data from file
exp_name = 'C:\Users\UserName\Documents\data\ExperimentData.lif';


current_dir = cd;
cd('C:\Users\UserName\Documents\MATLAB\bfmatlab');
data = bfopen(exp_name);
cd(current_dir)


%%
for i = 1:1:size(data,1)
    adress =  data{i,1}{1,2};
    ind_semicol = strfind(data{i,1}{1,2},';');
    name = adress(ind_semicol(1)+2:ind_semicol(2)-1);
    %name = adress(ind_semicol(1)+2:end);
    fprintf(['%d - ' name '\n'],i);
end

%%
Start_sequence_res = 2; %last line on res of last series
N_series = 1;%1:size(data,1) % for all the cells

int_range_fluo = [0 4];

adress =  data{N_series,1}{1,2};
ind_semicol = strfind(data{N_series,1}{1,2},';');
name = adress(ind_semicol(1)+2:ind_semicol(2)-1);


fprintf('--- Series # %d ---\n',N_series)

start_from = 25; %25 excude first 6 frames for fast bleaching

% load time series in matrix
T_tot = size(data{N_series,1},1)/2;
Im = data{N_series,1}{start_from,1};
DataMat_big = nan([size(Im),T_tot-start_from+1]);
DataMat_bf = nan([size(Im),T_tot-start_from+1]);

for i = start_from:T_tot
    Im = data{N_series,1}{i*2-1,1}; %fluorescence
    Im_bf = data{N_series,1}{i*2,1}; %fluorescence
    DataMat_big(:,:,i-start_from+1) = Im;
    DataMat_bf(:,:,i-start_from+1) = Im_bf;
end

%%

intTot = squeeze(nanmean(nanmean(DataMat_big,1),2));
bleach = (intTot(1)-intTot(end))/intTot(1)*100;
fprintf('Intensity first to last frame: -%.1f%% \n',bleach)

figure
IntTbf = squeeze(mean(mean(DataMat_bf,2),1));
%mean(IntTbf)
IntTbfScaled = (IntTbf-mean(IntTbf(1:20)))/mean(IntTbf(1:20))*100;
plot(IntTbfScaled)
xlabel('frame')
ylabel('laser intensity (%)')


Im_ave_big = mean(DataMat_big,3);
Im_ave_bf = mean(DataMat_bf,3);

% intensity range
[vals,bins_e] = hist(double(Im_ave_big(:)),100);
csum = cumsum(vals);
csum = csum./csum(end);
[~,I] = min(abs(csum-1));
Int_max = bins_e(I);



k = true;
while k
    
    figure
    imshow(Im_ave_bf,[],'initialmagnification',300)
    figure
    imshow(Im_ave_big,int_range_fluo,'initialmagnification',300)
    
    if size(res,1)>Start_sequence_res
        label = 1;
        for i = Start_sequence_res:size(res,1)
            pos = res{i,13};
            %rectangle('Position',pos, 'LineWidth',2,'EdgeColor','white')
            text(floor(pos(1)+pos(3)/3),floor(pos(2)+pos(3)/2),num2str(label),'FontSize',18,'Color','white');
            label = label+1;
        end
    end
    
    % get rectangle around single cell
    rect_cell = round(getrect());
    rect_cell(3:4) = min(rect_cell(3:4));
    DataMat_c = DataMat_big(rect_cell(2):rect_cell(2)+rect_cell(4),rect_cell(1):rect_cell(1)+rect_cell(3),:);
    DataMat_bf_c = DataMat_bf(rect_cell(2):rect_cell(2)+rect_cell(4),rect_cell(1):rect_cell(1)+rect_cell(3),:);
    close
    
    ShowNandBwithColormap(DataMat_c);
    
    pause(2)
    
    rotation = [0,0,0,0];
    DataMat = DataMat_c;
    
    %     % pause
    %     %
    % ------------------------------------------------------------------------%
        rotations = cell(1);
        rotations(1,1) = {[0,0,0,0]};
        rotations(1,2) = {[5,5,1,3]};
        rotations(1,3) = {[5,5,3,3]};
        rotations(1,4) = {[5,5,3,5]};
        rotations(1,5) = {[5,5,5,5]};
        rotations(1,6) = {[5,5,5,10]};
    
        out_compare = compare_cellmovement(DataMat_c,rotations);
    
        rotation = out_compare{2};
        DataMat = out_compare{1};
    % ------------------------------------------------------------------------%
    
    close all
    Im_ave_bf = mean(DataMat_bf,3);
    
    Im_ave = mean(DataMat,3);
    
    % intensity range
    [vals,bins] = hist(double(Im_ave(:)),100);
    csum = cumsum(vals);
    csum = csum./csum(end);
    
    [~,I] = min(abs(csum-1));
    Int_max = bins(I);
    
    B_s = [];
    N_s = [];
    e_s = [];
    n_s = [];
    e_box = [];
    n_box = [];
    rect_all5 = [];
    
    bleach5 = [];
    INT_bf5 = [];
    INT_fluo5 = [];
    
    bpall=[];
    
    for j=1:5 %average over 5 regions
        fprintf('box # %d \n',j)
        % get rectangle
        imshow(Im_ave,int_range_fluo,'initialmagnification',600)
        rect = round(getrect());
        close
        
        
        
        % matrix of selection
        DataMat_c = DataMat(rect(2):rect(2)+rect(4)+1,rect(1):rect(1)+rect(3)+1,:);
        
        % display time series of selection with slider
        %figure_slider(DataMat_c,Int_max)
        
        % Calculate numbers and brightness
        v = var(DataMat_c,0,3);
        m = mean(DataMat_c,3);
        
        
        B_p = v./m;
        N_p = m.^2./v;
        e_p = B_p-1;
        %n_p = (m.^2)./(v - m);
        n_p = N_p.*B_p./e_p;
        
        bpall=[bpall;e_p(:)];
        
        B_s = [B_s,median(B_p(:))];
        N_s = [N_s,median(N_p(:))];
        e_s = [e_s,median(e_p(:))];
        n_s = [n_s,median(n_p(:))];
        
        
        
        % N&B Boxcar        
        vBox = movvar(DataMat_c,50,0,3,'Endpoints','discard');
        mBox = movmean(DataMat_c,50,3,'Endpoints','discard');       
        B_pBox = nanmean(vBox./mBox,3);
        N_pBox = nanmean(mBox.^2./vBox,3);
        e_pBox = B_pBox-1;
        n_pBox = N_pBox.*B_pBox./e_pBox;
        
        
        e_box = [e_box,median(e_pBox(:))];
        n_box = [n_box,median(n_pBox(:))];
        rect_all5 = [rect_all5;rect];
        
        
        % for bleaching
        DataMat_bf_cc = DataMat_bf_c(rect(2):rect(2)+rect(4)+1,rect(1):rect(1)+rect(3)+1,:);
        intTot_bf = squeeze(nanmean(nanmean(DataMat_bf_cc,1),2));
        
        intTot_fluo = squeeze(nanmean(nanmean(DataMat_c,1),2));
        start_fluo = mean(intTot_fluo(1:10));
        end_fluo = mean(intTot_fluo(end-10:end));
        bleach1 = (start_fluo-end_fluo)/start_fluo*100;
        fprintf('Intensity first to last frame: -%.1f%% \n',bleach1)
        
        INT_bf5 = [INT_bf5,intTot_bf];
        INT_fluo5 = [INT_fluo5,intTot_fluo];
        
        bleach5 = [bleach5,bleach1];
    end
    
    %     hist(bpall(:),50)
    %     pause
    
    B = mean(B_s);
    B_std = std(B_s);
    N = mean(N_s);
    N_std = std(N_s);
    e = mean(e_s);
    e_std = std(e_s);
    n = mean(n_s);
    n_std = std(n_s);
    e_box = mean(e_box);
    n_box = mean(n_box);
    
    fprintf('Brightness (sigma/k): %2.4f +/- %2.4f \n',e,e_std);
    fprintf('Numbers (k^2/sigma): %2.4f +/- %2.4f \n',n,n_std);
    
    
    % pause %ok?
    
    % addline on storage
    
    d = size(res,1)+1;
    res(d,1) = {exp_name};
    res(d,2) = {'comment'}; % type of sample
    res(d,3) = {name};
    res(d,4) = {N};
    res(d,5) = {B};
    res(d,6) = {n};
    res(d,7) = {e};
    res(d,8) = {N_std};
    res(d,9) = {B_std};
    res(d,10) = {n_std};
    res(d,11) = {e_std};
    res(d,12) = {mean(bleach5)};
    res(d,13) = {rect_cell};
    res(d,14) = {rect_all5};
    res(d,15) = {rotation};
    res(d,16) = {INT_fluo5};
    res(d,17) = {INT_bf5};
    res(d,18) = {n_box};
    res(d,19) = {e_box};
end

