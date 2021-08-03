%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

clear;
close all;
% file with model
respath='./';
outpath='./Results/';

% number of periods and burn-in
NT_sim=5000;
NT_ini=round(NT_sim/10);

% EE error
compEEErr=1;

% Force the creation of a sim_res file
force_output = 1;
resfile_list={'res_20201117_base'}; 

% resfile_struct = dir('res_20201117*.mat');
% resfile_cell = {resfile_struct(:).name};
% ff=1;
% for f=1:length(resfile_cell)
%     resfile=resfile_cell{f};   
%     resfile_tokens =  strsplit(resfile,'.');
%     resfile_tokens = strsplit(resfile_tokens{1},'_');
%     resfile_new = [resfile_tokens{1},'_',resfile_tokens{2},'_',resfile_tokens{3}];
% %    if ~isfile(['sim_',resfile_new,'.mat'])
%     if length(resfile_tokens)>3
%         movefile(resfile, [resfile_new,'.mat']);
%         resfile_list{ff}=resfile_new;
%         ff=ff+1;
%     end
% end



for f=1:length(resfile_list)

resfile=resfile_list{f};   
load([respath,resfile]);


% output table file
outstats=[outpath,'stats_',resfile,'fit.xls'];
errstats=[outpath,'errstats_',resfile,'fit.xls'];

%% list of variable names to be reported (must match names defined in
% compStSt and in this file)

%--------------------------------------------------------------------------
% start simulation
%--------------------------------------------------------------------------

% set starting point
start_ex=5;
startpt=struct;
startpt.K=stv.State.K;
startpt.bC=stv.State.bC;
if ~mobj.Params.C_bank_only
    startpt.bS=stv.State.bS;
    startpt.KSsh=stv.State.KSsh;
end
startpt=orderfields(startpt,mobj.En_names);
startpt_vec=model.DSGEModel.structToVec(startpt)';
startpt_vec=[start_ex,startpt_vec];

% simulate
[simseries,varnames,errmat]=mobj.simulate(NT_sim,NT_ini,startpt_vec,compEEErr,[]);
simseries_orig=simseries;
varnames_orig=varnames;
statevec = simseries(:,1);

[simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames);
nvars = length(varnames);

% Create table object for easier access
simtable=array2table(simseries);
[~,ia,~]=unique(varnames);
simtable=simtable(:,ia);
simtable.Properties.VariableNames=varnames(ia);

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:nvars
    indexmap.put(varnames{i},i);
end



%% regression for epsilon calibration

if ~mobj.Params.C_bank_only
qS=simseries(:,indexmap.get('qS'));
qC=simseries(:,indexmap.get('qC'));
AC=simseries(:,indexmap.get('AC'));
AS=simseries(:,indexmap.get('AS'));
FS=simseries(:,indexmap.get('FS'));
FC=simseries(:,indexmap.get('FC'));
rf=simseries(:,indexmap.get('rf'));
varrho=simseries(:,indexmap.get('vartheta'));
GDP=simseries(:,indexmap.get('GDP'));
C=simseries(:,indexmap.get('C'));
SDf=simseries(:,indexmap.get('MRS_C'));
valtemp = unique(simtable.varrho);
idxR = simtable.varrho ==max(valtemp);
xt = simtable.LS(idxR).*simtable.varrho(idxR)./...
     simtable.lS(idxR);


MRS_C=simseries(:,indexmap.get('MRS_C'));
MRS_S=simseries(:,indexmap.get('MRS_S'));
FSzero=(FS==0);
logCY = log(C./GDP);

y=qC-qS;
%y=MRS_C-MRS_S;
%   X=[ log(AS./GDP) log(AC./GDP)  FS*100 logCY [0; diff(log(GDP))] ];
 X=[  log(AS./GDP) log(AC./GDP)    FS*100 logCY   [0; diff(log(GDP))]];
%X=[  log(AS./GDP) FS*100  ];
[EstCov,se,coeff]  = hac(X ,y);
AScoeff=coeff(2);
ACcoeff=coeff(3);
disp(['Coefficient on AS: ', num2str(AScoeff)]);
disp(['Coefficient on AC: ', num2str(ACcoeff)]);
disp(['SE on AS: ', num2str(se(2))]);
disp(['SE on AC: ', num2str(se(3))]);


disp(' ');

% add to simseries as constant
simseries=[simseries, ones(size(simseries,1),1)*AScoeff];
indexmap.put('AScoeff',size(simseries,2));
nvars=nvars+1;

simseries=[simseries, ones(size(simseries,1),1)*ACcoeff];
indexmap.put('ACcoeff',size(simseries,2));
nvars=nvars+1;

simseries=[simseries, ones(size(simseries,1),1)*100*mean(1-xt)];
indexmap.put('Haircut',size(simseries,2));
nvars=nvars+1;


simseries=[simseries, AS./(AS+AC)];
indexmap.put('ASsh',size(simseries,2));
nvars=nvars+1;

% [~, ACAS_cycle] = hpfilter( AC./AS, 1600);
simseries=[simseries,AC./AS];
indexmap.put('ACAS',size(simseries,2));
nvars=nvars+1;


varnames=[varnames,{'AScoeff', 'ACcoeff', 'Haircut','ASsh','ACAS'}];

end

%% --------------------------------------------------------------------------
% calculate stats
%--------------------------------------------------------------------------
varst=zeros(length(startpt_vec)-1,1);
for i=1:length(startpt_vec)-1
    varst(i)=indexmap.get(mobj.En_names{i});
end
% state variable means in stationary distribution
stvstat=mean(simseries(:,varst));

% calculate business cycle stats
% first for all periods, then separately for booms and recessions
smpsel={true(NT_sim-1,1), simseries(:,1) >= mobj.Params.muY, simseries(:,1) < mobj.Params.muY, simseries(:,3) > 0};
statsout=cell(numel(smpsel),1);

GDPind=indexmap.get('GDP');
for j=1:numel(smpsel)
    % subsample
    simtmp=simseries(smpsel{j},:);
    statstmp=zeros(nvars,11);
    statstmp(:,1)=nanmean(simtmp)';
    statstmp(:,2)=nanstd(simtmp)';
    % contemp and first-order autocorrelations
    autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
    conm=autocorrm(1:nvars,1:nvars);
    lagm=autocorrm(nvars+1:end,1:nvars);
    % corr with shocks
    statstmp(:,3:4)=[conm(:,1),lagm(1,:)'];
    statstmp(:,5:6)=[conm(:,2),lagm(2,:)'];
    statstmp(:,7:8)=[conm(:,3),lagm(3,:)'];
    statstmp(:,9:10)=[conm(:,GDPind),lagm(GDPind,:)'];
   % vector with fo autocorr
    statstmp(:,11)=diag(lagm);
    statsout{j}=statstmp;
end

%--------------------------------------------------------------------------
% output
%--------------------------------------------------------------------------

% overview output for eyeball check against analytic st.st. values
% make one big structure with steady-state values
stvbig=model.HelperCollection.combineStructs({stv.Sol,stv.State,stv.Add,stv.statsout});

% output table
% make index vector
dispnames=unique(varnames);
[displist,dispnames]=model.HelperCollection.makeListFromNames(indexmap,dispnames);
ndvars=length(displist);

disp(' ');
disp('Simulation steady state');

% overview output 
fprintf('Frequency: ');
for j=1:numel(smpsel)
    % select vars
    tabout{j}=statsout{j}(displist,:);
    fprintf('%f\t',sum(smpsel{j}));
end
fprintf('\n');
disp('-------------');

for s=1:ndvars
    if isfield(stvbig,dispnames{s})
        ststval=stvbig.(dispnames{s});
    else
        ststval=0;
    end
    if numel(dispnames{s}) > 7
        fprintf('%d\t%4s\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    else
        fprintf('%d\t%4s\t\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    end
    %fprintf('%d\t%4s\t\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    for j=1:numel(smpsel)
        fprintf('\t%f, %f |',tabout{j}(s,1),tabout{j}(s,2));
    end
    fprintf('\n');
end

% check grid bounds
if ~mobj.Params.C_bank_only
    state_range=5:8;
else
    state_range=5:6;
end
min_vec=min(simseries(:,state_range));
max_vec=max(simseries(:,state_range));
disp('State bounds:');
disp(mobj.Pfct.SSGrid.StateBounds(:,2:end));
disp('Simulation mins:');
disp(min_vec);
disp('Simulation max:');
disp(max_vec);

if compEEErr
    avg_err=mean(abs(errmat))';
    med_err=median(abs(errmat))';
    p75_err=prctile(abs(errmat),75)';
    p95_err=prctile(abs(errmat),95)';
    p99_err=prctile(abs(errmat),99)';
    p995_err=prctile(abs(errmat),99.5)';
    max_err=max(abs(errmat))';
    errtab=table(avg_err,med_err,p75_err,p95_err,p99_err,p995_err,max_err);
    errarr=table2array(errtab);
    disp(' ');
    disp('-----------------------------------------------');
    disp('Average and maximum Euler equation error');
    fprintf('Equ.no.\t\tAvg.\t\tMed.\t\tp75\t\t\tp95\t\t\tp99\t\t\tp99.5\t\tMax.\n');
    for s=1:length(avg_err)
        fprintf('%d\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,errarr(s,1),errarr(s,2),errarr(s,3), ...
            errarr(s,4),errarr(s,5),errarr(s,6),errarr(s,7));
    end
    
    % plot EE error for these equations
%     plotEE_pol=[3,4];
%     plotEE_state=[0,0];
%     for i=1:length(plotEE_pol)
%         points=simseries(:,[3,4]);
%         errvals=abs(errmat(1:end-1,plotEE_pol(i)));
%         if plotEE_state(i)>0
%             itmp=(statvec==plotEE_state(i));
%             points=points(itmp,:);
%             errvals=errvals(itmp,:);
%         end
%         model.HelperCollection.scatterPoints2D(points,errvals);
%     end
    
end

%%
% write to file
if usejava('desktop') 
    values=struct2cell(mobj.Params);
    paramout=cell2table(values,'RowNames',fieldnames(mobj.Params));
    colnames={'mean','std','corrY','corrY_1','corrZ','corrZ_1','corrVR','corrVR_1','corrGDP','corrGDP_1','AC'};
      for j=1:numel(tabout)
        tableout=array2table(tabout{j},'RowNames',dispnames,'VariableNames',colnames);
        writetable(tableout,outstats,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);
      end
    writetable(paramout,outstats,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');
    writetable(errtab,errstats,'FileType','spreadsheet');
end    

params=mobj.Params;
disp(['Saving simulation data to .mat file: ',['sim_',resfile,'.mat']]);
save(['sim_',resfile,'.mat'],'simseries','displist','dispnames','errmat','tabout','outstats', ...
    'errstats','errtab','indexmap','NT_ini','NT_sim','smpsel','statevec','statsout','varnames','params');

clearvars={'simseries','displist','dispnames','errmat','tabout','outstats', ...
    'errstats','errtab','indexmap','smpsel','statevec','statsout','varnames','params'};
clear(clearvars{:});

end
