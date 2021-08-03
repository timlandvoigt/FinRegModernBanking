if usejava('desktop')
   clear; 
end
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
outpath='./Results/';
if ~exist('resfile','var')
    resfile_list={'res_20201117_base'};
end
% for transition dynamics load starting point from base case
resfile_startpt=[];
%resfile_startpt='res_20200605_base';
no_par_processes = 16;

for f=1:length(resfile_list)
    
    resfile=resfile_list{f};    
    
    load([respath,resfile,'.mat']);
    
    % Initial Economy Config
    varlist={'simseries','statevec','indexmap','varnames'};
    if ~isempty(resfile_startpt)
        % just simulate unconditional dynamics
        load(['sim_',resfile_startpt],varlist{:});
        %start_shock=0;
        start_shock=[0,1,2];
    else        
        % run shock comparison
        load(['sim_',resfile],varlist{:});
        start_shock=[0,1,2];
    end

    % number of periods and burn-in
    N_shock=length(start_shock);
    N_runs=10000;
    NT_sim=20;
    NT_ini=0;
    NT_sim=NT_sim+1;    
    
    % set starting point
    start_ini=9;
    startpt_mat=[];
    if ~isempty(resfile_startpt)
        statevec=statevec(2:end);
        simfrac=histc(statevec,1:mobj.Exogenv.exnpt)/length(statevec);
        for s=1:length(simfrac)
            startpt=struct;
            startvals=mean(simseries(statevec==s,:));
            startvals=startvals(1:end-4); % no regression coefs
            N_vars=length(startvals);
            startpt.K=startvals(indexmap.get('K'));
            startpt.bS=startvals(indexmap.get('bS'));
            startpt.bC=startvals(indexmap.get('bC'));
            startpt.KSsh=startvals(indexmap.get('KSsh'));
            startpt=orderfields(startpt,mobj.En_names);
            startpt_vec=model.DSGEModel.structToVec(startpt)';
            startpt_vec=[s,startpt_vec];           
            startpt_mat=[startpt_mat; repmat(startpt_vec,floor(N_runs*simfrac(s)),1)];
        end
    end
    if size(startpt_mat,1)<N_runs
        startpt=struct;
        statevec=statevec(2:end);
        startvals=mean(simseries(statevec==start_ini,:));
        startvals=startvals(1:end-5);% no regression coefs
        N_vars=length(startvals);
        startpt.K=startvals(indexmap.get('K'));
        startpt.bS=startvals(indexmap.get('bS'));
        startpt.bC=startvals(indexmap.get('bC'));
        startpt.KSsh=startvals(indexmap.get('KSsh'));
        startpt=orderfields(startpt,mobj.En_names);
        startpt_vec=model.DSGEModel.structToVec(startpt)';
        startpt_vec=[start_ini,startpt_vec];
        startpt_mat=[startpt_mat; repmat(startpt_vec,N_runs-size(startpt_mat,1),1)];
    end
        
    % compute Euler equation error?
    compEEErr=1;
    % make graphs grayscale
    grayscale=0;
    
    % output table file
    outfile=['GTR_',resfile];
    
    simseries_median = cell(N_shock,1);
    simseries_mean = cell(N_shock,1);
    simseries_std = cell(N_shock,1);
    
    simseries_diff_median = cell(N_shock,1);
    simseries_diff_mean = cell(N_shock,1);
    simseries_diff_std = cell(N_shock,1);
    
    open_parpool;
    
    
    for s=1:N_shock
        
        disp(['Shock ',num2str(s),' of ',num2str(N_shock)]);
        tens_simseries = zeros(NT_sim,N_vars,N_runs);
        
        % compute entry of random number matrix that sets first state
        % deterministically to start_shock
        if start_shock(s)>0
            transprob=cumsum(mobj.Exogenv.mtrans(start_ini,:));
            shock_prob=transprob(start_shock(s));
            if start_shock(s)>1
                shock_prob_minus=transprob(start_shock(s)-1);
            else
                shock_prob_minus=0;
            end
            rvar_next=(shock_prob+shock_prob_minus)/2;
        end
        
        % Create shock matrix
        rng(150,'combRecursive');
        %shmatfull = rand(NT_sim*N_runs,1);
        shmatfull = rand(N_runs,NT_sim);
        
        fprintf([repmat('.',1,100) '\n\n']);
        
        parfor n=1:N_runs
%        for n=1:N_runs            
            %--------------------------------------------------------------------------
            % start simulation
            %--------------------------------------------------------------------------
            %fprintf('Run %d - Start \n',n);
            % simulate
            shmat = shmatfull(n,:)';
            if start_shock(s)>0
                shmat(1)=rvar_next;
            end
            startvec=startpt_mat(n,:);
            [simseries,varnames]=mobj.simulate(NT_sim,NT_ini,startvec,compEEErr,shmat);
            simseries_orig=simseries;
            varnames_orig=varnames;
            statevec = simseries(:,1);
            %fprintf('Run %d - After simulation \n',n);
            
            [simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames);
            nvars = length(varnames);
            %fprintf('Run %d - After computation \n',n);
            
            tens_simseries(:,:,n) = [startvals; simseries];
            
            if mod(n,N_runs/100)==0
                %disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
                fprintf('\b|\n');
            end
            
        end
        fprintf('\n');

        
        simseries_median{s} = median(tens_simseries,3);
        simseries_mean{s} = mean(tens_simseries,3);
        simseries_std{s} = std(tens_simseries,[],3);
        
        if start_shock(s)==0 || start_shock(s)==start_ini
            % If no shock, store for later differencing
            tens_simseries_0 = tens_simseries;
        else
            % If actual shock, difference and save
            tens_simseries_diff = tens_simseries - tens_simseries_0;
            simseries_diff_median{s} = median(tens_simseries_diff,3);
            simseries_diff_mean{s} = mean(tens_simseries_diff,3);
            simseries_diff_std{s} = std(tens_simseries_diff,[],3);
        end
        
    end
    
        
    save(outfile,'simseries_mean','simseries_median','simseries_std', ...
        'simseries_diff_mean', 'simseries_diff_median', 'simseries_diff_std', ...
        'indexmap','NT_sim','N_shock');    
    
end
