clear,clc
%% read data
[data_g_se,name]=xlsread('cag abundance.xlsx',1);
id_genus=name(1,2:end);
genus_se=name(2:end,1);

seledata=sum(data_g_se,2);
x=find(seledata>0);
data_g_se=data_g_se(x,:);
genus_se=genus_se(x);

%% sample id
[~,id_aexm0]=xlsread('sample annotation.xlsx',8);
[c,ia,~]=intersect(id_genus,id_aexm0);
id=c;
datase=data_g_se(:,ia)';

%%
[R,P]=corr(datase,'type','spearman');
padj=zeros(size(P));
for i=1:size(P,1)
    padj(i,:)=mafdr(P(i,:));
end

%%
[a,b]=find(P<0.05);
pselect=zeros(size(a,1),2);
for i=1:size(a,1)
if a(i)<b(i) 
    pselect(i,:)=[a(i) b(i)];
else pselect(i,:)=[0 0];
end
end
pselect(all(pselect==0,2),:)=[];

rselect=zeros(size(pselect,1),5);
for j=1:size(pselect,1)
    if abs(R(pselect(j,1),pselect(j,2)))>0
      rselect(j,:)=[pselect(j,1),pselect(j,2),R(pselect(j,1),pselect(j,2)),P(pselect(j,1),pselect(j,2)),padj(pselect(j,1),pselect(j,2))];  
    else rselect(j,:)=[0 0 0 0 0];
    end
end
rselect(all(rselect==0,2),:)=[];

out_var=[genus_se(rselect(:,1)) genus_se(rselect(:,2))];
out_coefficient=rselect(:,[3 4 5]);
%%
auni=unique(a);