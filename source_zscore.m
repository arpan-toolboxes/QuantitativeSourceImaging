function [Z_score_c,Z_score_p] =source_zscore(source_dics_stim,model,z_thr)
%%

source_P=source_dics_stim.pow((source_dics_stim.inside),1);
source_c=source_dics_stim.pos(source_dics_stim.inside,:);
Temp=[source_P source_c];

[values, order] = sort(Temp(:,2));
SortedTemp = Temp(order,:);

mid=[];
for ii=1:size(SortedTemp,1)
    if (SortedTemp(ii,2)==0 && SortedTemp(ii,3)==0 && SortedTemp(ii,4)==0) ;
        mid=ii;
    end
end

LB=SortedTemp(1:mid,:);
RB=SortedTemp(mid+1:size(SortedTemp,1),:);

LB_source_P=LB(:,1);
LB_source_c=LB(:,2:4);

RB_source_P=RB(:,1);
RB_source_c=RB(:,2:4);

Center=[];
j=1;
for aa=1:size(LB_source_c,1)
    if abs(LB_source_c(aa,1))<=15
        Center(j)=aa;
        j=j+1;
    end
end
LB_source_c(Center,:)=[];
LB_source_P(Center,:)=[];

Center=[];
j=1;
for aa=1:size(RB_source_c,1)
    if abs(RB_source_c(aa,1))<=15
        Center(j)=aa;
        j=j+1;
    end
end
RB_source_c(Center,:)=[];
RB_source_P(Center,:)=[];

Std=std(LB_source_P);
Mean=mean(LB_source_P);
for jj=1:length(LB_source_P)
    Zscore_DICS_LB(jj)=(LB_source_P(jj)-Mean)/Std;
end

Std=std(RB_source_P);
Mean=mean(RB_source_P);
for jj=1:length(RB_source_P)
    Zscore_DICS_RB(jj)=(RB_source_P(jj)-Mean)/Std;
end

for zz=1:length(z_thr)
    Z_DICS_LB=[];
    Z_DICS_LB_P=[];
    Z_stat=quantile(Zscore_DICS_LB,z_thr(zz));
    Z_ind=find(Zscore_DICS_LB>Z_stat);
    Z_DICS_LB=LB_source_c(Z_ind,:);
    Z_DICS_LB_P=LB_source_P(Z_ind,:);
    
    Z_DICS_RB=[];
    Z_DICS_RB_P=[];
    Z_stat=quantile(Zscore_DICS_RB,z_thr(zz));
    Z_ind=find(Zscore_DICS_RB>Z_stat);
    Z_DICS_RB=RB_source_c(Z_ind,:);
    Z_DICS_RB_P=RB_source_P(Z_ind,:);
    
    Max_Vox=[];
    if  strcmp(model,'dist')==1 || strcmp(model,'two')==1  
        Z_score_c{zz}=[Z_DICS_LB; Z_DICS_RB];
        Z_score_p{zz}=[Z_DICS_LB_P; Z_DICS_RB_P];
    elseif strcmp(model,'single')==1
        Z_score_c{zz}=[Z_DICS_LB];
        Z_score_p{zz}=[Z_DICS_LB_P];
        
    end
end