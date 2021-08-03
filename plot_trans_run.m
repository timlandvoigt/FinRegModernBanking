clear;
close all;

respath='./';
outpath='./Results/';
prefix='GTR';
plotshock=3;

if ~exist('econ1', 'var'), econ1='20201117_base'; end
if ~exist('econ2', 'var'), econ2='20201117_nu10'; end
if ~exist('colors', 'var'), colors = {'k-o', 'b-o'}; end
%if ~exist('label', 'var'), label = {'\theta=17%', '\theta=25%'}; end
if ~exist('label', 'var'), label = {}; end

% Check number
if exist('econ3', 'var')
    N_economy = 3;
    econlist={econ1,econ2,econ3};
elseif exist('econ2', 'var')
    N_economy=2; % numbers economies to plot on same graph
    econlist={econ1,econ2};
else
    N_economy=1; % numbers economies to plot on same graph
    econlist={econ1};
end

outfile=[prefix,'_res_run_'];

simseries_diffmean_econ=cell(N_economy,1);
simseries_mean_econ=cell(N_economy,1);
for n=1:N_economy
    outfile=[outfile,econlist{n}];
    resfile=['res_',econlist{n}];
    load([respath,resfile,'.mat']);
    load([respath,'sim_',resfile,'.mat']);
    load([respath,prefix,'_',resfile,'.mat']);
        
    simseries_diffmean_econ{n} = simseries_diff_mean{plotshock};
    simseries_mean_econ{n} =  simseries_mean{plotshock};
    clear simseries_mean;   
end

   
tvec=0:NT_sim-1;   

%% Make other Graphs

% file names for graphs (set to empty for no printing)
printfiles={[outpath,outfile,'IRF1'],[outpath,outfile,'IRF2'],[outpath,outfile,'IRF3'],[outpath,outfile,'IRF4']};        

% which variables
brsel1=[indexmap.get('Z'),indexmap.get('I'),indexmap.get('GDP'),...
        indexmap.get('C'),indexmap.get('H'),indexmap.get('DWL')];  
brsel2=[indexmap.get('KSsh'),indexmap.get('Slev'),indexmap.get('Clev'),...
        indexmap.get('FS'),indexmap.get('FC'),indexmap.get('rateS')];  

brsel_all=[brsel1,brsel2];

relative_plot=[ones(1,5),2,ones(1,6)];


% transforms
plain = @(x)(x);
prct = @(x)(100*x);
prct_change = @(x)((x/x(1) -1)*100);
transformvec={plain,plain,plain,plain,plain,plain,...
              plain,plain,plain,plain,plain,plain};


nvar=length(brsel_all);   

brseries_gr=zeros(N_economy, NT_sim, nvar);
for s=1:N_economy
    for v=1:nvar
        thisv=brsel_all(v); 
        if relative_plot(v)==1
            this_series = simseries_diffmean_econ{s}(1:NT_sim,thisv)./simseries_mean_econ{s}(1:NT_sim,thisv) * 100;
        elseif relative_plot(v)==2
            this_series = simseries_diffmean_econ{s}(1:NT_sim,thisv) * 100;
        else
            this_series = simseries_mean_econ{s}(1:NT_sim,thisv);
        end
        this_series = transformvec{v}(this_series);
        brseries_gr(s,:,v) = this_series;
    end
end


titles1={'Z', 'BDS Investment', 'BDS Output', 'Consumption',  'Liquidity', 'DWL/GDP'}; %,['Discretionary spending ',lsuff]};
titles2={'S Capital Share', 'S Leverage', 'C Leverage', 'S Default Rate', 'C Default Rate', 'S Deposit Rate'}; %,['Discretionary spending ',lsuff]};

if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:6),tvec,titles1,colors,[2,3],label,printfiles{1});   
    makeImpRes(brseries_gr(:,:,7:12),tvec,titles2,colors,[2,3],label,printfiles{2});  
    %makeImpRes(brseries_gr(:,:,13:end),tvec,titles3,colors,[2,4],label,printfiles{3});  
end

