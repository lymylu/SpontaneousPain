BaseDir='S1optostimulus'; 
filelist=dir(BaseDir);
Eventtype={'PL5Hz','PL10Hz','PL20Hz','PLmix','S15Hz','S110Hz','S120Hz','S1mix'};
for i=3:length(filelist)
      tmpmat=matfile(fullfile(BaseDir,filelist(i).name),'Writable',true);
      tmpmat=cal_binspikes(tmpmat,0.1,Eventtype);
     
    try
      tmpmat=cal_optotag(tmpmat);
    end
      %tmpmat=cal_response(tmpmat,Eventtype);
    try
      tmpmat=cal_durationresponse(tmpmat,'Spontaneouspain','Spontaneousnopain');
    end
% tmpmat=cal_spectrogram(tmpmat);
end

%% summarized the data
clear PL* S1* Opto*
fileindex=[3,5,6];
varnames={'PL10Hz','PL5Hz','PL20Hz','Optotag','S110Hz','S15Hz','S120Hz','Spontaneouspain'};
for i=1:3
    tmpmat=matfile(fullfile(BaseDir,filelist(fileindex(i)).name));
    for j=1:length(varnames)
        eval([varnames{j},'(i)=NeuroResult(tmpmat.',varnames{j},');']);
    end
end
S15Hz_all=S15Hz.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','binspike'},[2,2,2,2,2,2]);
S110Hz_all=S110Hz.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','binspike'},[2,2,2,2,2,2]);
S120Hz_all=S120Hz.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','binspike'},[2,2,2,2,2,2]);
Optoresponse_all=Optotag.CollectSpikeVariables({'Optoresponse'},2);
Spontaneous_all=Spontaneouspain.CollectSpikeVariables({'response'},2);
Optoresponse=Optoresponse_all.Optoresponse;
Painresponse=cell2mat(Spontaneous_all.response);
valid=S110Hz_all.SPKinfo.firingRate>0.1&S110Hz_all.SPKinfo.troughToPeak>0.5;
tabulate(Optoresponse(valid));
tabulate(Painresponse(valid));

lightoffbefore=cell2mat(S110Hz_all.binspike);



function tmpdata=cal_durationresponse(tmpdata,yesdataname,nodataname)
    % cal the spike responses to the Spontaneous behavior
    % the classification of responses using the comparision of the true
    % firingrate of goal segements among that of the surrogate segments.
    nodata=eval(['tmpdata.',nodataname,';']);
    yesdata=eval(['tmpdata.',yesdataname,';']);
%     frequency=sum(cellfun(@(x) length(x),yesdata.SPKdata,'UniformOutput',1))/sum(yesdata.timerange(:,2)-yesdata.timerange(:,1));
%     nofrequency=sum(cellfun(@(x) length(x),nodata.spiketime,'UniformOutput',1))/sum(nodata.werange(:,2)-nodata.timerange(:,1));
    % permutation test?
    % nopain data and pain data
    nodata=NeuroResult(nodata);
    yesdata=NeuroResult(yesdata);
    nodata=nodata.Splice2Split;
    yesdata=yesdata.Splice2Split;
    nodata=nodata.NeuroResult2Struct;
    yesdata=yesdata.NeuroResult2Struct;
    notimerange=[min(reshape(nodata.SPKinfo.timerange,[],1)),max(reshape(nodata.SPKinfo.timerange,[],1))];
    yestimerange=[min(reshape(yesdata.SPKinfo.timerange,[],1)),max(reshape(yesdata.SPKinfo.timerange,[],1))];
    timerangeall=[min(min(notimerange),min(yestimerange)),max(max(notimerange),max(yestimerange))];
    yestimeall=sum(yesdata.SPKinfo.timerange(:,2)-yesdata.SPKinfo.timerange(:,1)); % the total time of spontaneouspain 
    notimeall=timerangeall(2)-timerangeall(1)-yestimeall; % the total time of spontaneousnopain
    notime=cell(1,size(nodata.SPKdata,1));
    yestime=cell(1,size(yesdata.SPKdata,1));
    for i=1:size(nodata.SPKdata,1)
        for j=1:size(nodata.SPKdata,2)
            notime{i}=cat(1,notime{i},nodata.SPKdata{i,j});
        end
        nofrequency{i}=length(notime{i})/notimeall;
    end
    for i=1:size(yesdata.SPKdata,1)
        for j=1:size(yesdata.SPKdata,2)
            yestime{i}=cat(1,yestime{i},yesdata.SPKdata{i,j});
        end
        yesfrequency{i}=length(yestime{i})/yestimeall;
    end
    spiketime=cellfun(@(x,y) sort(cat(1,x,y)),yestime,notime,'UniformOutput',0); % the spike time list from the whole time duration
    permutationpaintimerange=getPermutationtimerange(yesdata.SPKinfo.timerange,nodata.SPKinfo.timerange,notimeall,500);
    [response,permutationfreq]=cellfun(@(x,y) cal_permutationresponse(x,permutationpaintimerange,yestimeall,y),spiketime,yesfrequency,'UniformOutput',0);
    nodata.durationfrequency=nofrequency;
    yesdata.durationfrequency=yesfrequency;
    yesdata.response=response;
    yesdata.permutationfreq=permutationfreq;
    nodata.timeall=spiketime;
    eval(['tmpdata.',nodataname,'=nodata;']);
    eval(['tmpdata.',yesdataname,'=yesdata;']);
end
function datamat=cal_binspikes(datamat,Fs,fieldname) 

for c=1:length(fieldname)
tmpmat=eval(['datamat.',fieldname{c},';']);
try
spiketime=arrayfun(@(x) x{:}+tmpmat.EVTinfo.timerange(1),tmpmat.SPKdata,'UniformOutput',0);
catch
    spiketime=arrayfun(@(x) x{:}-tmpmat.EVTinfo.timestart(1),tmpmat.SPKdata,'UniformOutput',0);
end
binspike=[];
for k=1:size(spiketime,1)
    for i=1:size(spiketime,2)
    st.spiketime=spiketime{k,i};
    try
    [binspike{k}(:,i),t]=binspikes(st.spiketime,Fs,tmpmat.EVTinfo.timerange);
    catch
         [binspike{k}(:,i),t]=binspikes(st.spiketime,Fs,[0,round(tmpmat.EVTinfo.timestop-tmpmat.EVTinfo.timestart)]);
    end
    end
end
    tmpmat.binspike=binspike;
    tmpmat.binspiket=t;
    eval(['datamat.',fieldname{c},'=tmpmat;']);
end
end% using the psth from chronux toolbox ().
function permutationpaintimerange=getPermutationtimerange(yestimerange,notimerange,notimeall,nperm)
        rng("default");
        for i=1:nperm
        seqofyes=randperm(size(yestimerange,1)); % randperm the sequence of pain time segment
%         notimeseq=diff([sort(randsample(1:notimeall,length(yestimerange)+1)),notimeall(end)]); % randperm the time sequences between permutation spontanouespain.
        notimeseq=diff(sort([notimeall*rand(size(yestimerange,1)+2,1)]));
        yesduration=yestimerange(seqofyes,2)-yestimerange(seqofyes,1); % each time duration of the permutation spontaneouspain.
        timerangeall=[min(min(notimerange),min(yestimerange)),max(max(notimerange),max(yestimerange))];
        for j=1:length(notimeseq)
            if j==2
            permutationpaintimerange(j-1,1,i)=sum(notimeseq(1:j-1));
            permutationpaintimerange(j-1,2,i)=permutationpaintimerange(j-1,1,i)+yesduration(j-1);
            elseif j>2
                permutationpaintimerange(j-1,1,i)=sum(notimeseq(1:j-1))+sum(yesduration(1:j-2));
                permutationpaintimerange(j-1,2,i)=permutationpaintimerange(j-1,1,i)+yesduration(j-1);
            end
        end
        permutationpaintimerange(:,:,i)=permutationpaintimerange(:,:,i)+timerangeall(1);
        if permutationpaintimerange(j-1,2,i)>timerangeall(end)
            a=1;
        end
    end
end
function [response,permutationfreq]=cal_permutationresponse(spiketime,permutationpaintimerange,yesduration,yesfrequency)
        permutationfreq=[];
        for i=1:size(permutationpaintimerange,3) 
            permutationspike=[];
        for j=1:size(permutationpaintimerange,1)
            permutationspike=cat(1,permutationspike,length(find(spiketime>permutationpaintimerange(j,1,i)&spiketime<permutationpaintimerange(j,2,i))));
        end
            permutationfreq(i)=sum(permutationspike)/sum(yesduration);
        end
    p=length(find(permutationfreq>yesfrequency))/size(permutationpaintimerange,3);
    if ~isempty(spiketime)
        if p>=0.95
        response=-1;
        elseif p<=0.05
        response=1;
        else
        response=0;
        end
    else
        response=0;
    end
end
function datamat=cal_optotag(datamat)
    tmpmat=datamat.Optotag;
    index=ismember(tmpmat.EVTinfo.eventdescription,'OptoTag');
    tmpspiketime=arrayfun(@(x) x{:}+tmpmat.EVTinfo.timerange(1),tmpmat.SPKdata(:,index),'UniformOutput',0);
    for k=1:size(tmpspiketime)
        for i=1:size(tmpspiketime,2)
            st.spiketime=tmpspiketime{k,i};
            tmpbinspike{k}(:,i)=binspikes(st.spiketime,1000,tmpmat.EVTinfo.timerange);
        end
    end
    for i=1:length(tmpbinspike)
    t=linspace(-2,4,size(tmpbinspike{i},1));    
    p(i)=salt(tmpbinspike{i}(t<0&t>-0.5,:)',tmpbinspike{i}(t>0,:)',0.001,0.005);
    if p(i)<0.05
        tmpmat.Optoresponse(i)=1;
    else
        firingbefore=mean(tmpbinspike{i}(t<0&t>-0.5,:),1);
        firingafter=mean(tmpbinspike{i}(t>0&t<0.5,:),1);
        [~,p_tmp]=ttest(firingbefore,firingafter);
        if p_tmp<0.05 && mean(firingbefore)>mean(firingafter)
            tmpmat.Optoresponse(i)=-1;
        else
            tmpmat.Optoresponse(i)=0;
        end
    end
    end
    datamat.Optotag=tmpmat;
end