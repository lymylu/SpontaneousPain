%% Spontaneous Analysis of Spike
active='Spontaneouspain';
noactive='Spontanousnopain';
basepath='Spikefiles/';
Group={'(\<RAT6)|(\<RAT11)|(\<RAT14)|(\<RAT16)','(\<RAT4)|(\<RAT7)|(\<RAT8)|(\<RAT15)|(\<RAT17)|(\<RAT19)|(\<RAT21)'};
groupname={'Saline','CFA'};
Time={'(_1d_)|(_3d_)','(_5d_)|(_7d_)|(_9d_)|(_10d_)|(_11d_)|(_14d_)'}; % baseline ,early stage late stage
Brainregion={'(_cluster2_)|(_cluster1_)|(_cluster3_)}(_cluster4_)|(_cluster5_)|(_cluster6_)|(_cluster7_)|(_cluster8_)'};% IL cluster1,2 PL cluster4-6
brainname={'PL'};
spikelist=dir(basepath);
spikelist=struct2table(spikelist);
spikelist=spikelist.name;
% Pyramidal neuron
PL=spikeselect(basepath,spikelist,active,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% PL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{2});
resultfieldname=fieldnames(PL);
for i=1:length(resultfieldname)
     tmp1=eval(['PL.',resultfieldname{i}]);
     tmp2=[];
%      tmp2=eval(['PL2.',resultfieldname{i}]);
   
    try
        eval(['SpontaneousPain.PL_PY_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
    catch
        eval(['SpontaneousPain.PL_PY_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
    end
end
% IL=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% IL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{1});
% resultfieldname=fieldnames(IL);
% for i=1:length(resultfieldname)
%     tmp1=eval(['IL.',resultfieldname{i}]);
%     tmp2=eval(['IL2.',resultfieldname{i}]);
% %     tmp2=[];
%     try
%         eval(['SpontaneousPain.IL_PY_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
%     catch
%         eval(['SpontaneousPain.IL_PY_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
%     end
% end
%% Interneuron
PL=spikeselect(basepath,spikelist,active,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% PL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{2});
resultfieldname=fieldnames(PL);
for i=1:length(resultfieldname)
    tmp1=eval(['PL.',resultfieldname{i}]);
%     tmp2=eval(['PL2.',resultfieldname{i}]);
     tmp2=[];
    try
        eval(['SpontaneousPain.PL_IN_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
    catch
        eval(['SpontaneousPain.PL_IN_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
    end
end
%[IL_nopain_In,ILname_nopain_In,ILbasefiringrate_In]=spikeselect(basepath,spikelist,nopain,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{2});
% IL=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% IL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{1});
% resultfieldname=fieldnames(IL);
% for i=1:length(resultfieldname)
%     tmp1=eval(['IL.',resultfieldname{i}]);
%     tmp2=eval(['IL2.',resultfieldname{i}]);
% %     tmp=[];
%    try
%         eval(['SpontaneousPain.IL_IN_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
%     catch
%         eval(['SpontaneousPain.IL_IN_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
%     end
% end
%%
figure; dataresponse={'SpontaneousPain.PL_PY_response','SpontaneousPain.IL_PY_response','SpontaneousPain.PL_IN_response','SpontaneousPain.IL_IN_response'};
basefiring={'SpontaneousPain.PL_PY_firingrate','SpontaneousPain.IL_PY_firingrate','SpontaneousPain.PL_IN_firingrate','SpontaneousPain.IL_IN_firingrate'};
titleresponse={'PL pyramidal','IL pyramidal','PL interneuron','IL interneuron'};
for i=1:2:3
subplot(2,2,i)
tmp=eval(dataresponse{i});
tmpbase=eval(basefiring{i});
% basename=strrep(dataresponse{i},'response','');
% basename2=strrep(basename,'SpontaneousPain','');
% [index]=findsub(eval(['SpontaneousPain',basename2,'spikename']),eval(['Movement',basename2,'spikename']),eval(['find(Movement',basename2,'firingrate>0.1)']));
% tmp=tmp(index);
%      tmp(tmpbase<0.1)=[];
% %      tmpbase(tmpbase<0.1)=[];
%      firingResponse{i}{1}=tmpbase(find(tmp==1));
%      firingResponse{i}{2}=tmpbase(find(tmp==-1));
%      firingResponse{i}{3}=ttest2(firingResponse{i}{1},firingResponse{i}{2});
x=tabulate(tmp); 
xin=x(find(x(:,1)==-1),2);
xex=x(find(x(:,1)==1),2);
xno=x(find(x(:,1)==0),2);
if isempty(xex) xex=0;end;
if isempty(xno) xno=0;end;
if isempty(xin) xin=0;end;
pie([xin,xno,xex],{['inhibitory N=',num2str(xin)],['no change N=',num2str(xno)],['excitatory N=',num2str(xex)]});
title(titleresponse{i});
end
PYfiring={'SpontaneousPain.PL_PY_nopainfiringrate','SpontaneousPain.PL_PY_painfiringrate','SpontaneousPain.IL_PY_nopainfiringrate','SpontaneousPain.IL_PY_painfiringrate'};
INfiring={'SpontaneousPain.PL_IN_nopainfiringrate','SpontaneousPain.PL_IN_painfiringrate','SpontaneousPain.IL_IN_nopainfiringrate','SpontaneousPain.IL_IN_painfiringrate'};
basefiring={'SpontaneousPain.PL_PY_firingrate','SpontaneousPain.PL_PY_firingrate','SpontaneousPain.IL_PY_firingrate','SpontaneousPain.IL_PY_firingrate'};
for i=1:2
    tmp=eval(PYfiring{i});
 tmpbase=eval(basefiring{i});
      tmp(tmpbase<0.1)=[];
    PYfiringrate{i}=tmp;
end
figure;
subplot(1,2,1);
bar([mean(PYfiringrate{1}),mean(PYfiringrate{2})]);ylim([0,3]);
 [~,p1]=ttest(PYfiringrate{1},PYfiringrate{2});
% [~,p2]=ttest(PYfiringrate{3},PYfiringrate{4});
xticklabels({'PLnopain','PLpain'}); title('Pyramidal neuron');
basefiring={'SpontaneousPain.PL_IN_firingrate','SpontaneousPain.PL_IN_firingrate','SpontaneousPain.IL_IN_firingrate','SpontaneousPain.IL_IN_firingrate'};
for i=1:2
    tmp=eval(INfiring{i});
 tmpbase=eval(basefiring{i});
      tmp(tmpbase<0.1)=[];
    INfiringrate{i}=tmp;
end
subplot(1,2,2); bar([mean(INfiringrate{1}),mean(INfiringrate{2})]);ylim([0,4]);
[~,p3]=ttest(INfiringrate{1},INfiringrate{2});
% [~,p4]=ttest(INfiringrate{3},INfiringrate{4});
xticklabels({'PLnopain','PLpain'});title('Interneuron');

%%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
% % % %  move and nomove comparison
active='move';
nomove='nomove';
basepath='LFP_spike/';
Group={'(\<RAT25)|(\<RAT26)'}; % cfa RAT4,7,8, saline RAT11,14,6
groupname={'CFA','Saline'};
Time={'(_5d_)|(_7d_)|(_9d_)|(_10d_)|(_11d_)|(_14d_)'}; % baseline ,early stage late stage
Brainregion={'(_cluster2_)|(_cluster1_)|(_cluster3_)}(_cluster4_)|(_cluster5_)|(_cluster6_)|(_cluster7_)|(_cluster8_)'};% IL cluster1,2 PL cluster4-6
brainname={'IL','PL'};
spikelist=dir(basepath);
spikelist=struct2table(spikelist);
spikelist=spikelist.name;
% Pyramidal neuron
PL=spikeselect(basepath,spikelist,active,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% PL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{2});
resultfieldname=fieldnames(PL);
for i=1:length(resultfieldname)
     tmp1=eval(['PL.',resultfieldname{i}]);
     tmp2=[];
%      tmp2=eval(['PL2.',resultfieldname{i}]);
   
    try
        eval(['Movement.PL_PY_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
    catch
        eval(['Movement.PL_PY_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
    end
end
% IL=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% % IL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{1});
% resultfieldname=fieldnames(IL);
% for i=1:length(resultfieldname)
%     tmp1=eval(['IL.',resultfieldname{i}]);
% %     tmp2=eval(['IL2.',resultfieldname{i}]);
%      tmp2=[];
%     try
%         eval(['Movement.IL_PY_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
%     catch
%         eval(['Movement.IL_PY_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
%     end
% end
% Interneuron
PL=spikeselect(basepath,spikelist,active,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
%  PL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{2});
resultfieldname=fieldnames(PL);
for i=1:length(resultfieldname)
    tmp1=eval(['PL.',resultfieldname{i}]);
%     tmp2=eval(['PL2.',resultfieldname{i}]);
     tmp2=[];
    try
        eval(['Movement.PL_IN_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
    catch
        eval(['Movement.PL_IN_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
    end
end
%[IL_nopain_In,ILname_nopain_In,ILbasefiringrate_In]=spikeselect(basepath,spikelist,nopain,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{2});
% IL=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1});
% IL2=spikeselect(basepath,spikelist,move,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{1});
% resultfieldname=fieldnames(IL);
% for i=1:length(resultfieldname)
%     tmp1=eval(['IL.',resultfieldname{i}]);
%     tmp2=eval(['IL2.',resultfieldname{i}]);
% %     tmp=[];
%    try
%         eval(['Movement.IL_IN_',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
%     catch
%         eval(['Movement.IL_IN_',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
%     end
% end
%%
figure; dataresponse={'Movement.PL_PY_response','Movement.IL_PY_response','Movement.PL_IN_response','Movement.IL_IN_response'};
basefiring={'Movement.PL_PY_firingrate','Movement.IL_PY_firingrate','Movement.PL_IN_firingrate','Movement.IL_IN_firingrate'};
titleresponse={'PL pyramidal','IL pyramidal','PL interneuron','IL interneuron'};
for i=1:2:3
subplot(2,2,i)
tmp=eval(dataresponse{i});
tmpbase=eval(basefiring{i});
    tmp(tmpbase<0.1)=[];
x=tabulate(tmp); 
xin=x(find(x(:,1)==-1),2);
xex=x(find(x(:,1)==1),2);
xno=x(find(x(:,1)==0),2);
if isempty(xex) xex=0;end;
if isempty(xno) xno=0;end;
if isempty(xin) xin=0;end;
pie([xin,xno,xex],{['inhibitory N=',num2str(xin)],['no change N=',num2str(xno)],['excitatory N=',num2str(xex)]});
title(titleresponse{i});
end
PYfiring={'Movement.PL_PY_nomovefiringrate','Movement.PL_PY_movefiringrate','Movement.IL_PY_nomovefiringrate','Movement.IL_PY_movefiringrate'};
INfiring={'Movement.PL_IN_nomovefiringrate','Movement.PL_IN_movefiringrate','Movement.IL_IN_nomovefiringrate','Movement.IL_IN_movefiringrate'};
basefiring={'Movement.PL_PY_firingrate','Movement.PL_PY_firingrate','Movement.IL_PY_firingrate','Movement.IL_PY_firingrate'};
for i=1:2
    tmp=eval(PYfiring{i});
    tmpbase=eval(basefiring{i});
      tmp(tmpbase<0.1)=[];
    PYfiringrate{i}=tmp;
end
figure;
subplot(1,2,1);
bar([mean(PYfiringrate{1}),mean(PYfiringrate{2})]);ylim([0,3]);
[~,p1]=ttest(PYfiringrate{1},PYfiringrate{2});
% [~,p2]=ttest(PYfiringrate{3},PYfiringrate{4});
xticklabels({'PLnomove','PLmove'}); title('Pyramidal neuron');
basefiring={'Movement.PL_IN_firingrate','Movement.PL_IN_firingrate','Movement.IL_IN_firingrate','Movement.IL_IN_firingrate'};
for i=1:2
    tmp=eval(INfiring{i});
       tmpbase=eval(basefiring{i});
      tmp(tmpbase<0.1)=[];
    INfiringrate{i}=tmp;
end
subplot(1,2,2); bar([mean(INfiringrate{1}),mean(INfiringrate{2})]);ylim([0,4]);
[~,p3]=ttest(INfiringrate{1},INfiringrate{2});
% [~,p4]=ttest(INfiringrate{3},INfiringrate{4});
xticklabels({'PLnomove','PLmove'});title('Interneuron');
%%
% % % % 
% response venn between spontanoues response and move response
[index1,index2]=findsub(SpontaneousPain.PL_PY_spikename,Movement.PL_PY_spikename,Movement.PL_PY_firingrate>0.1);
SpontaneousresponsePL=SpontaneousPain.PL_PY_response(index1);
MovementresponsePL=Movement.PL_PY_response(index2);
x=find(SpontaneousresponsePL==1&MovementresponsePL==1);
z=find(SpontaneousresponsePL==-1&MovementresponsePL==-1);
% [index1,index2]=findsub(SpontaneousPain.IL_PY_spikename,Movement.IL_PY_spikename,Movement.IL_PY_firingrate>0.1);
% SpontaneousresponseIL=SpontaneousPain.IL_PY_response(index1);
% MovementresponseIL=Movement.IL_PY_response(index2);
% y=find(SpontaneousresponseIL==-1&MovementresponseIL==-1);
% w=find(SpontaneousresponseIL==1&MovementresponseIL==1);
figure;
subplot(2,2,1);
venn([length(find(SpontaneousresponsePL==1)),length(find(MovementresponsePL==1))],length(x));
disp(num2str([length(find(SpontaneousresponsePL==1)),length(find(MovementresponsePL==1)),length(x)]));
title('the overlay of PL excitatory responses in SpontaneousPain and Movement');
% subplot(2,2,2);
% venn([length(find(SpontaneousresponseIL==-1)),length(find(MovementresponseIL==-1))],length(y));
% disp(num2str([length(find(SpontaneousresponseIL==-1)),length(find(MovementresponseIL==-1)),length(y)]));
% title('the overlay of IL inhibitory responses in SpontaneousPain and Movement');
subplot(2,2,3);
venn([length(find(SpontaneousresponsePL==-1)),length(find(MovementresponsePL==-1))],length(z));
disp(num2str([length(find(SpontaneousresponsePL==-1)),length(find(MovementresponsePL==-1)),length(z)]));
title('the overlay of PL inhibitory responses in SpontaneousPain and Movement');
% subplot(2,2,4);
% venn([length(find(SpontaneousresponseIL==1)),length(find(MovementresponseIL==1))],length(w));
% disp(num2str([length(find(SpontaneousresponseIL==1)),length(find(MovementresponseIL==1)),length(w)]));
% title('the overlay of IL excitatory responses in SpontaneousPain and Movement');

%%
% % response venn between laser evoked response and spontanoues response;
PL_PY_Laser=spikeselectcombined({basepath,spikelist,'VonfreyLow','firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1}},[]);
% IL_PY_Laser=spikeselectcombined({basepath,spikelist,'Laser','firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1}},...
% {basepath,spikelist,'VonfreyLow','firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{1}});

PL_IN_Laser=spikeselectcombined({basepath,spikelist,'VonfreyLow','firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1}},[]);
% IL_IN_Laser=spikeselectcombined({basepath,spikelist,'Laser','firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time{1},'Brainregion',Brainregion{1}},...
% {basepath,spikelist,'Laser','firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time{1},'Brainregion',Brainregion2{1}});
%%
[index1,index2]=findsub(PL_PY_Laser.spikename,SpontaneousPain.PL_PY_spikename,SpontaneousPain.PL_PY_firingrate>0.1);
laserresponsePL_PY=PL_PY_Laser.response(index1);
spontaneousPL_PY=SpontaneousPain.PL_PY_response(index2);
% [index1,index2]=findsub(IL_PY_Laser.spikename,SpontaneousPain.IL_PY_spikename,SpontaneousPain.IL_PY_firingrate>0.1);
% laserresponseIL_PY=IL_PY_Laser.response(index1);
% spontaneousIL_PY=SpontaneousPain.IL_PY_response(index2);
[index1,index2]=findsub(PL_IN_Laser.spikename,SpontaneousPain.PL_IN_spikename,SpontaneousPain.PL_IN_firingrate>0.1);
laserresponsePL_IN=PL_IN_Laser.response(index1);
spontaneousPL_IN=SpontaneousPain.PL_IN_response(index2);
% [index1,index2]=findsub(IL_IN_Laser.spikename,SpontaneousPain.IL_IN_spikename,SpontaneousPain.IL_IN_firingrate>0.1);
% laserresponseIL_IN=IL_IN_Laser.response(index1);
% spontaneousIL_IN=SpontaneousPain.IL_IN_response(index2);
%  x1=find(laserresponseIL_PY==1&spontaneousIL_PY==1);
 y=find(laserresponsePL_PY==-1&spontaneousPL_PY==-1);
%  z=find(laserresponseIL_PY==-1&spontaneousIL_PY==-1);
 w=find(laserresponsePL_PY==1&spontaneousPL_PY==1);
 figure;
 subplot(2,2,1);
 x=tabulate(laserresponsePL_PY); 
    xin=x(find(x(:,1)==-1),2);
    xex=x(find(x(:,1)==1),2);
    xno=x(find(x(:,1)==0),2);
    pie([xin,xno,xex],{['inhibitory N=',num2str(xin)],['no change N=',num2str(xno)],['excitatory N=',num2str(xex)]});
%  subplot(2,2,2);
%   x=tabulate(laserresponseIL_PY); 
%     xin=x(find(x(:,1)==-1),2);
%     xex=x(find(x(:,1)==1),2);
%     xno=x(find(x(:,1)==0),2);
%     pie([xin,xno,xex],{['inhibitory N=',num2str(xin)],['no change N=',num2str(xno)],['excitatory N=',num2str(xex)]});
subplot(2,2,3)
    x=tabulate(laserresponsePL_IN); 
    xin=x(find(x(:,1)==-1),2);
    xex=x(find(x(:,1)==1),2);
    xno=x(find(x(:,1)==0),2);
    pie([xin,xno,xex],{['inhibitory N=',num2str(xin)],['no change N=',num2str(xno)],['excitatory N=',num2str(xex)]});
%  subplot(2,2,4);
%   x=tabulate(laserresponseIL_IN); 
%     xin=x(find(x(:,1)==-1),2);
%     xex=x(find(x(:,1)==1),2);
%     xno=x(find(x(:,1)==0),2);
%     pie([xin,xno,xex],{['inhibitory N=',num2str(xin)],['no change N=',num2str(xno)],['excitatory N=',num2str(xex)]});    

figure;

subplot(2,2,1);
venn([length(find(laserresponsePL_PY==1)),length(find(spontaneousPL_PY==1))],length(w));
disp(num2str([length(find(laserresponsePL_PY==1)),length(find(spontaneousPL_PY==1)),length(w)]));
title('the overlay of PL excitatory responses in SpontaneousPain and Laser');
% subplot(2,2,2);
% venn([length(find(laserresponseIL_PY==-1)),length(find(spontaneousIL_PY==-1))],length(z));
% disp(num2str([length(find(laserresponseIL_PY==-1)),length(find(spontaneousIL_PY==-1)),length(z)]));
% title('the overlay of IL inhibitory responses in SpontaneousPain and Laser');
subplot(2,2,3);
venn([length(find(laserresponsePL_PY==-1)),length(find(spontaneousPL_PY==-1))],length(y));
disp(num2str([length(find(laserresponsePL_PY==-1)),length(find(spontaneousPL_PY==-1)),length(y)]));
title('the overlay of PL inhibitory responses in SpontaneousPain and Laser');
% subplot(2,2,4);
% venn([length(find(laserresponseIL_PY==1)),length(find(spontaneousIL_PY==1))],length(x1));
% disp(num2str([length(find(laserresponseIL_PY==1)),length(find(spontaneousIL_PY==1)),length(x1)]));
% title('the overlay of IL excitatory responses in SpontaneousPain and Laser');

 
%%
% show the typical response neuron inhibitory IL pain #86
% excitatory PL pain # 98
Spontaneousshow(SpontaneousPain.PL_PY_spikename{98},SpontaneousPain.PL_PY_painfiringrate(98),'Spontanousnopain','Spontaneouspain');
 xlim([800,1300]);
Spontaneousshow(SpontaneousPain.IL_PY_spikename{91},SpontaneousPain.IL_PY_painfiringrate(91),'Spontanousnopain','Spontaneouspain');
xlim([200,900]);
%% cal the basefiring rate in different time
Time2={'(_baseline_)|(_BASELINE_)|(_Baseline_)','(_1d_)','(_3d_)','(_7d_)|(_9d_)','(_14d_)|(_15d_)'}; 
for i=1:length(Time2)
    PL=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time2{i},'Brainregion',Brainregion{2});
    PL2=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time2{i},'Brainregion',Brainregion2{2});
    PL_PY_basefiring{i}=cat(2,PL.firingrate,PL2.firingrate);
    IL=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group{1},'Time',Time2{i},'Brainregion',Brainregion{1});
    IL2=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Pyramidal Cell','Group',Group2{1},'Time',Time2{i},'Brainregion',Brainregion2{1});
    IL_PY_basefiring{i}=cat(2,IL.firingrate,IL2.firingrate);
    PL=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time2{i},'Brainregion',Brainregion{2});
    PL2=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time2{i},'Brainregion',Brainregion2{2});
    PL_IN_basefiring{i}=cat(2,PL.firingrate,PL2.firingrate);
    IL=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group{1},'Time',Time2{i},'Brainregion',Brainregion{1});
    IL2=basefiringselect(basepath,spikelist,'firingrate',[0,Inf],'celltype','Narrow Interneuron','Group',Group2{1},'Time',Time2{i},'Brainregion',Brainregion2{1});
    IL_IN_basefiring{i}=cat(2,IL.firingrate,IL2 .firingrate);
end
PL_PY_base=cellfun(@(x) x(x>0.5),PL_PY_basefiring,'UniformOutput',0);
IL_PY_base=cellfun(@(x) x(x>0.5),IL_PY_basefiring,'UniformOutput',0);


    function Spontaneousshow(spikename,spikefrequency,nofieldname,fieldname)
tmpmat=matfile(spikename);
[~,permutationfreq,spiketime]=Spontaneous_response_cal(tmpmat,spikefrequency,nofieldname,fieldname);
figure;
subplot(1,2,1);
a=hist(permutationfreq,20);
hist(permutationfreq,20);
hold on;
line([spikefrequency,spikefrequency],[0,max(a)]);
axis tight;
ylabel('Permutation counts');
xlabel('Spike frequency');
subplot(1,2,2);
binsize=10;
binrange=linspace(min(spiketime),max(spiketime),(max(spiketime)-min(spiketime))/binsize+1);
a=hist(spiketime,binrange);
hist(spiketime,binrange);
tmp=tmpmat.Spontaneouspain;
tmptimerange=tmp.timerange;
hold on
for i=1:size(tmptimerange,1)
    fill([tmptimerange(i,1),tmptimerange(i,1),tmptimerange(i,2),tmptimerange(i,2)],[0,max(a),max(a),0],'r','FaceAlpha',.3,'EdgeAlpha',0);
end
title(spikename);
line([min(spiketime),max(spiketime)],[mean(permutationfreq)*binsize,mean(permutationfreq)*binsize]);
ylabel('Spike counts per bin');
xlabel('Time (s)'); axis tight;

end
    function [Coherence,t,f]=Coherence_cal(LFP1,LFP2)
        params.pad=0;
        params.Fs=1250;
        params.fpass=[0,100];
        params.err=0;
        params.trialave=1;
        params.tapers=[3,5];
        params.segave=1;
        [Coherence,~,~,~,t,f]=coherencysegc(squeeze(mean(LFP1,2)),squeeze(mean(LFP2,2)),2,params);
    end
    function [pdcatb,pdcbta,freq]=PDC_cal(LFP1,LFP2)
    % cal the Granger Causuality of LFP1 and LFP2 using the GCCA toolbox
        LFP1=squeeze(mean(LFP1,2));
        LFP2=squeeze(mean(LFP2,2));
        data=cat(2,LFP1,LFP2);
        data=downsample(data,5);
        timeindex=windowepoched(data,[2,2],0,2,250);
        for i=1:size(timeindex,2)
            dataepoch(:,:,i)=data(timeindex(:,i),:);
        end
        dataepoch=permute(dataepoch,[2,1,3]);
        [aic,bic,moaic,mobic] = tsdata_to_infocrit(dataepoch,20,'LWR',false);
       [A,SIG,E] = tsdata_to_var(dataepoch,moaic,'LWR');
       [G,info] = var_to_autocov(A,SIG);
       [f,fres] = autocov_to_spwcgc(G,256);
        pdcatb=squeeze(f(2,1,:));
        pdcbta=squeeze(f(1,2,:));
        pdcatb=mean(pdcatb,2);
        pdcbta=mean(pdcbta,2);
        freq=sfreqs(256,250);
    end 
    function Result=spikeselect(varargin)
% usage: the optional parameters is firingrate, celltype(Pyramidal Cell,
% Narrow Interneuron or Wide Interneuron), Group (string compare of rat
% num), Brain region(string compare of clusters), and Time (string compare
% of time);


p=inputParser;
addRequired(p,'basepath');
addRequired(p,'spikelist');
addRequired(p,'eventtype');
addParameter(p,'firingrate',[0,Inf],@isnumeric);
addParameter(p,'celltype','',@(x) ischar(x));
addParameter(p,'Group',[],@(x) ischar(x));
addParameter(p,'Brainregion',[],@(x) ischar(x));
addParameter(p,'Time',[],@(x) ischar(x));
parse(p,varargin{:});
if strcmp(p.Results.eventtype,'Spontaneouspain')
    Result.paintimeall=[];Result.nopaintimeall=[];Result.painfiringrate=[];Result.nopainfiringrate=[];Result.response=[];
end
if strcmp(p.Results.eventtype,'move')
    Result.movetimeall=[];Result.nomovetimeall=[];Result.movefiringrate=[];Result.nomovefiringrate=[];Result.response=[];
end
spikelist=p.Results.spikelist;
Filefilter={'Group','Brainregion','Time'};
k=1;
    for i=1:length(Filefilter)
        tmp=eval(['p.Results.',Filefilter{i},';']);
        if ~isempty(tmp)
          index(:,k)=cellfun(@(x) ~isempty(regexpi(x,tmp,'match')),spikelist,'UniformOutput',1);
        else
          index(:,k)=ones(length(spikelist),1);
        end
        k=k+1;
    end
    index=prod(index,2);
    spike=spikelist(logical(index));
    n=1;
    multiWaitbar(['Load eventtype:',p.Results.eventtype],0);
    for i=1:length(spike)
        try
        tmp=matfile(fullfile(p.Results.basepath,spike{i}));
        celltype=tmp.celltype;
        try
           celltype=tmp.celltype1;
        end
        if contains(celltype,p.Results.celltype)&&(tmp.Firingrate>p.Results.firingrate(1)&&tmp.Firingrate<p.Results.firingrate(2))
               eventdata=eval(['tmp.',p.Results.eventtype,';']);
               Result.frequency(n)=sum(cellfun(@(x) length(x),eventdata.spiketime,'UniformOutput',1))/sum(eventdata.timerange(:,2)-eventdata.timerange(:,1));
               if strcmp(p.Results.eventtype,'Spontaneouspain')
                  [Result.response(n),~,~,Result.paintimeall(n),Result.nopaintimeall(n),Result.nopainfiringrate(n),Result.painfiringrate(n)]=Spontaneous_response_cal(tmp,Result.frequency(n),'Spontanousnopain','Spontaneouspain');
               elseif strcmp(p.Results.eventtype,'move')
                   [Result.response(n),~,~,Result.movetimeall(n),Result.nomovetimeall(n),Result.nomovefiringrate(n),Result.movefiringrate(n)]=Spontaneous_response_cal(tmp,Result.frequency(n),'nomove','move');
               else
                   try
                    Result.response(n)=eventdata.response_2;
                   catch
                       Result.response(n)=eventdata.response;
                   end
               end
               Result.spikename{n}=tmp.Properties.Source;
               Result.firingrate(n)=tmp.Firingrate;
               n=n+1;
               try
                   Result.firingrate(n)=tmp.Firingrate1;
               catch
                   a=1;
               end
           
        end
        end
        multiWaitbar('LoadSpikes...',i/length(spike));
    end
    multiWaitbar('LoadSpikes','close');
    multiWaitbar(['Load eventtype:',p.Results.eventtype],1);
    end
    
    function [response,permutationfreq,spiketime,paintimeall,nopaintimeall,nopainfiringrate,painfiringrate]=Spontaneous_response_cal(tmpdata,frequency,nofieldname,fieldname)
    % cal the spike responses to the Spontaneous behavior
    nodata=eval(['tmpdata.',nofieldname]);
    yesdata=eval(['tmpdata.',fieldname]);
    timeall=sum(nodata.timerange(:,2)-nodata.timerange(:,1));
    spiketime=[];
    % permutation test criteria
    nopaintimerange=[min(reshape(nodata.timerange,[],1)),max(reshape(nodata.timerange,[],1))];
    paintimerange=[min(reshape(yesdata.timerange,[],1)),max(reshape(yesdata.timerange,[],1))];
    timerangeall=[min(min(nopaintimerange),min(paintimerange)),max(max(nopaintimerange),max(paintimerange))];
    paintimeall=sum(yesdata.timerange(:,2)-yesdata.timerange(:,1)); % the total time of spontaneouspain 
    nopaintimeall=timerangeall(2)-timerangeall(1)-paintimeall;% the total time of spontaneousnopain
    spiketime=[];
    for i=1:length(nodata.spiketime)
        spiketime=cat(1,spiketime,nodata.spiketime{i});
    end
    nopainfiringrate=length(spiketime)/nopaintimeall;
    spiketime1=[];
    for i=1:length(yesdata.spiketime)
        spiketime1=cat(1,spiketime1,yesdata.spiketime{i});
    end
    painfiringrate=length(spiketime1)/paintimeall;
    spiketime=sort(cat(1,spiketime,spiketime1)); % the spike time list from the whole time duration
    nperm=500;
    % resample the data from the combined nopain and pain time all;
%     for i=1:nperm
%         seqofpain=randperm(length(yesdata.spiketime)); % randperm the sequence of pain time segment
%         nopaintimeseq=diff([sort(randsample(1:nopaintimeall*1000,length(yesdata.spiketime)+1))./1000,nopaintimeall(end)]); % randperm the time sequences between permutation spontanouespain.
%         painduration=yesdata.timerange(seqofpain,2)-yesdata.timerange(seqofpain,1); % each time duration of the permutation spontaneouspain.
%         for j=1:length(nopaintimeseq)
%             if j>1
%             permutationpaintimerange(j-1,1)=sum(nopaintimeseq(1:j-1));
%             permutationpaintimerange(j-1,2)=permutationpaintimerange(j-1,1)+painduration(j-1);
%             elseif j>2
%                 permutationpaintimerange(j-1,1)=sum(nopaintimeseq(1:j-1))+sum(painduration(1:j-2));
%                 permutationpaintimerange(j-1,2)=permutationpaintimerange(j-1,1)+painduration(j-1);
%             end
%         end 
%         permutationpaintimerange=permutationpaintimerange+timerangeall(1);
%         permutationspike=[];
%         for j=1:length(permutationpaintimerange)
%             permutationspike=cat(1,permutationspike,length(find(spiketime>permutationpaintimerange(j,1)&spiketime<permutationpaintimerange(j,2))));
%         end
%         permutationfreq(i)=sum(permutationspike)/sum(painduration);
%     end
%     p=length(find(permutationfreq>frequency))/nperm;
%     if p>0.95
%         response=-1;
%     elseif p<0.05
%         response=1;
%     else
%         response=0;
%     end
%     only resample the nopain data.
    nopainspike=[];
    for i=1:size(nodata.timerange,1)
        if i==1
        nopainspike=cat(1,nopainspike,nodata.spiketime{i}-nodata.timerange(i,1));
        else
            timeperiod=nodata.timerange(:,2)-nodata.timerange(:,1);
            nopainspike=cat(1,nopainspike,nodata.spiketime{i}-nodata.timerange(i,1)+sum(timeperiod(1:i-1),1));
        end
    end
    nperm=500;
    for i=1:nperm
            begintime=randsample(0:(nopaintimeall-paintimeall)*1000,1)/1000;
            endtime=begintime+paintimeall;
            permutationfreq(i)=length(find(nopainspike>begintime&nopainspike<endtime))/paintimeall;  
    end
    p=length(find(permutationfreq>frequency))/nperm;
    if p>0.95
        response=-1;
    elseif p<0.05
        response=1;
    else
        response=0;
    end
    if (response==1 && nopainfiringrate>painfiringrate) || (response==-1 && nopainfiringrate<painfiringrate)
        response=0;
    end
    end
    function Result=basefiringselect(varargin)
    Result.firingrate=[];Result.spikename=[];
    p=inputParser;
    addRequired(p,'basepath');
    addRequired(p,'spikelist');
    addParameter(p,'firingrate',[0,Inf],@isnumeric);
    addParameter(p,'celltype','',@(x) ischar(x));
    addParameter(p,'Group',[],@(x) ischar(x));
    addParameter(p,'Brainregion',[],@(x) ischar(x));
    addParameter(p,'Time',[],@(x) ischar(x));
    parse(p,varargin{:});
    spikelist=p.Results.spikelist;
    Filefilter={'Group','Brainregion','Time'};
    k=1;
    for i=1:length(Filefilter)
        tmp=eval(['p.Results.',Filefilter{i},';']);
        if ~isempty(tmp)
          index(:,k)=cellfun(@(x) ~isempty(regexpi(x,tmp,'match')),spikelist,'UniformOutput',1);
        else
          index(:,k)=ones(length(spikelist),1);
        end
        k=k+1;
    end
    index=prod(index,2);
    spike=spikelist(logical(index));
    n=1;
    for i=1:length(spike)
        try
        tmp=matfile(fullfile(p.Results.basepath,spike{i}));
        celltype=tmp.celltype;
        try
           celltype=tmp.celltype1;
        end
        if contains(celltype,p.Results.celltype)&&(tmp.Firingrate>p.Results.firingrate(1)&&tmp.Firingrate<p.Results.firingrate(2))
               Result.spikename{n}=tmp.Properties.Source;
               Result.firingrate(n)=tmp.Firingrate;
               n=n+1;
        end
        end
    end
    end
    function Result=spikeselectcombined(input1,input2)
    Result=[];
    Tmp=spikeselect(input1{:}); 
    if ~isempty(input2)
        Tmp2=spikeselect(input2{:});
    end
    resultfieldname=fieldnames(Tmp);
    for i=1:length(resultfieldname)
        tmp1=eval(['Tmp.',resultfieldname{i}]);
        if ~isempty(input2)
            tmp2=eval(['Tmp2.',resultfieldname{i}]);
        else
            tmp2=[];
        end
        try
            eval(['Result.',resultfieldname{i},'=cat(1,tmp1,tmp2);']);
        catch
            eval(['Result.',resultfieldname{i},'=cat(2,tmp1,tmp2);']);
        end
    end
    end
    function [index1,index2]=findsub(A,B,restrict)
    % find the same spikename  in A from B
    % restrict is the condition of B
    for i=1:length(A)
        [~,A{i}]=fileparts(A{i});
    end
    
    for i=1:length(B)
        [~,B{i}]=fileparts(B{i});
    end
  C=B(restrict);
    for i=1:length(C)
        try
        index1(i)=find(strcmp(A,C{i})==1);
        catch
            index1(i)=nan;
        end
        index2(i)=find(strcmp(B,C{i})==1);
    end
    index2(isnan(index1))=[];
    index1(isnan(index1))=[];
    end