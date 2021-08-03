
clear; 
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
outpath='./Results/';
% for transition dynamics load starting point from base case
resfile_start='res_20201117_precrisis';
%resfile_end='res_20201117_postcrisis08';
resfile_end='res_20201117_postcrisis';
resfile_steps={'res_20201117_postcrisis085',...
               'res_20201117_postcrisis085',...
               'res_20201117_postcrisis09',...
               'res_20201117_postcrisis09',...
               'res_20201117_postcrisis095',...
               'res_20201117_postcrisis095',...
               'res_20201117_postcrisis10',...
               'res_20201117_postcrisis10',...
               'res_20201117_postcrisis105',...
               'res_20201117_postcrisis105'};
% resfile_steps={'res_20201117_postcrisisc1',...
%                'res_20201117_postcrisisc1',...
%                'res_20201117_postcrisisc2',...
%                'res_20201117_postcrisisc2',...
%                'res_20201117_postcrisisc3',...
%                'res_20201117_postcrisisc3',...
%                'res_20201117_postcrisisc4',...
%                'res_20201117_postcrisisc4',...
%                'res_20201117_postcrisisc5',...
%                'res_20201117_postcrisisc5'};
%resfile_steps=[];
% set starting point
% Initial Economy Config
varlist={'simseries','statevec','indexmap','varnames'};
% just simulate unconditional dynamics
load(['sim_',resfile_start],varlist{:});
mobjstart=load([respath,resfile_start,'.mat'],'mobj');
mobjstart=mobjstart.mobj;

% initial shock sequence before crash
start_shock_sequence=[2];
N_shock=length(start_shock_sequence);

% number of periods and burn-in
N_runs=5000;
NT_sim=30;
capdestruct=1;
no_par_processes=16;

start_ini=7;
startpt_mat=[];
statevec=statevec(2:end);
simfrac=histc(statevec,1:mobjstart.Exogenv.exnpt)/length(statevec);
for s=1:length(simfrac)
    startpt=struct;
    startvals=mean(simseries(statevec==s,:));
    startvals=startvals(1:end-4); % no regression coefs
    N_vars=length(startvals);
    startpt.K=startvals(indexmap.get('K'));
    startpt.bS=startvals(indexmap.get('bS'));
    startpt.bC=startvals(indexmap.get('bC'));
    startpt.KSsh=startvals(indexmap.get('KSsh'));
    startpt=orderfields(startpt,mobjstart.En_names);
    startpt_vec=model.DSGEModel.structToVec(startpt)';
    startpt_vec=[s,startpt_vec];
    startpt_mat=[startpt_mat; repmat(startpt_vec,floor(N_runs*simfrac(s)),1)];
end

if size(startpt_mat,1)<N_runs
    startpt=struct;
    statevec=statevec(2:end);
    startvals=mean(simseries(statevec==start_ini,:));
    startvals=startvals(1:end-4);% no regression coefs
    N_vars=length(startvals);
    startpt.K=startvals(indexmap.get('K'));
    startpt.bS=startvals(indexmap.get('bS'));
    startpt.bC=startvals(indexmap.get('bC'));
    startpt.KSsh=startvals(indexmap.get('KSsh'));
    startpt=orderfields(startpt,mobjstart.En_names);
    startpt_vec=model.DSGEModel.structToVec(startpt)';
    startpt_vec=[start_ini,startpt_vec];
    startpt_mat=[startpt_mat; repmat(startpt_vec,N_runs-size(startpt_mat,1),1)];
end

% save first row
N_varsim=mobjstart.NSTEX+mobjstart.NSTEN+mobjstart.NSOL+mobjstart.NV+mobjstart.NADD+10;
firstrow=[start_ini,startvals(1:N_varsim)];

% final periods
mobjfinal=load([respath,resfile_end,'.mat'],'mobj');
mobjfinal=mobjfinal.mobj;

% intermediate transition steps
N_steps=length(resfile_steps);
mobjlist=cell(N_steps,1);
for i=1:N_steps
    mobjthis=load([respath,resfile_steps{i},'.mat'],'mobj');
    mobjlist{i}=mobjthis.mobj;
end



% compute Euler equation error?
compEEErr=1;
% make graphs grayscale
grayscale=0;

% output table file
outfile=['GFC_',resfile_end];
   
%disp(['Shock ',num2str(s),' of ',num2str(N_shock)]);
tens_simseries = zeros(NT_sim,N_vars,N_runs);


% Create shock matrix
rng(150,'combRecursive');
shmatfull = rand(N_runs,NT_sim-N_shock);

open_parpool;

fprintf([repmat('.',1,100) '\n\n']);

parfor n=1:N_runs
%for n=1:N_runs
    %--------------------------------------------------------------------------
    % start simulation
    %--------------------------------------------------------------------------
    %fprintf('Run %d - Start \n',n);
    startvec=startpt_mat(n,:);
  
    % simulate
    shmat = shmatfull(n,:)';
    
    % compute entries of random number matrix that set first state
    % deterministically to start_shock_sequence
    rvar_sequence=zeros(N_shock,1);
    start_shock=startvec(1);
    for s=1:N_shock
        transprob=cumsum(mobjstart.Exogenv.mtrans(start_shock,:));
        shock_prob=transprob(start_shock_sequence(s));
        if start_shock_sequence(s)>1
            shock_prob_minus=transprob(start_shock_sequence(s)-1);
        else
            shock_prob_minus=0;
        end
        rvar_sequence(s)=(shock_prob+shock_prob_minus)/2;
        start_shock=start_shock_sequence(s);
    end
    
    
    % first period: run plus capital destruction
    [simseries_first,varnames,~,nextst]=mobjstart.simulate(N_shock,0,startvec,compEEErr,rvar_sequence,capdestruct);
    simseries_first_all = mobjstart.computeSimulationMoments(simseries_first,varnames,firstrow);
    firstperiod=simseries_first(end,:);
    N_varstotal=size(simseries_first_all,2);
    
    simseries_steps=zeros(N_steps,N_varstotal);
    for i=1:N_steps
        [simseries_this,varnames,~,nextst]=mobjlist{i}.simulate(1,0,nextst,compEEErr,shmat(i));   
        simseries_this_all = mobjlist{i}.computeSimulationMoments([firstperiod; simseries_this],varnames);
        simseries_steps(i,:)=simseries_this_all;   
        firstperiod=simseries_this(end,:);
    end
    
    [simseries_rest,varnames]=mobjfinal.simulate(NT_sim-N_shock-N_steps,0,nextst,compEEErr,shmat(N_steps+1:end));  
    [simseries_rest_all, varnames] = mobjfinal.computeSimulationMoments([firstperiod; simseries_rest],varnames);
    
%    tens_simseries(:,:,n) = [startvals; simseries_first_all; simseries_steps; simseries_rest_all];
    tens_simseries(:,:,n) = [simseries_first_all; simseries_steps; simseries_rest_all];
    
    if mod(n,N_runs/100)==0
        %disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
        fprintf('\b|\n');
    end
    
end
fprintf('\n');
%         varnames = varnames_store;
%         nvars = length(varnames);
%
%         % make HashMap with mapping of names to indices
%         indexmap=java.util.HashMap;
%         for i=1:nvars
%             indexmap.put(varnames{i},i);
%         end
%         varst=zeros(length(startpt_vec)-1,1);
%         for i=1:length(startpt_vec)-1
%             varst(i)=indexmap.get(mobj.En_names{i});
%         end

%save(outfile,'tens_simseries','indexmap');

simseries_median = median(tens_simseries,3);
simseries_mean = mean(tens_simseries,3);
simseries_std = std(tens_simseries,[],3);


save(outfile,'simseries_mean','simseries_median','simseries_std', ...
    'indexmap','NT_sim','N_shock');


