clear;
close all;

respath='./';
outpath='./Results/';
prefix='GFC';
datapath='Post Crisis.xlsx';

if ~exist('econ1', 'var'), econ1='20201117_postcrisis'; end
if ~exist('econ2', 'var'), econ2='20201117_postcrisis08'; end
if ~exist('colors', 'var'), colors = {'k--','k-o', 'b:'}; end
if ~exist('label', 'var'), label = {'data','\theta^{Post}=11%', '\theta^{Post}=8%'}; end

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

simseries_mean_econ=cell(N_economy,1);
for n=1:N_economy
    outfile=[outfile,econlist{n}];
    resfile=['res_',econlist{n}];
    load([respath,resfile,'.mat']);
    load([respath,'sim_',resfile,'.mat']);
    load([respath,prefix,'_',resfile,'.mat']);
        

    simseries_mean_econ{n} =  simseries_mean;

    clear simseries_mean;   
end

   
tvec=0:NT_sim-1;   

%% Make other Graphs

% file names for graphs (set to empty for no printing)
printfiles={[outpath,outfile,'IRF1'],[outpath,outfile,'IRF2'],[outpath,outfile,'IRF3'],[outpath,outfile,'IRF4']};        

% which variables
brsel4=[indexmap.get('GDP'),indexmap.get('C'),indexmap.get('H'),...
        indexmap.get('BSsh'),indexmap.get('Slev'),indexmap.get('Clev')];
    
brsel_all=[brsel4];

% transforms
plain = @(x)(x);
prct = @(x)(100*x);
prct_change = @(x)((x/x(1) -1)*100);
transformvec={prct_change,prct_change,prct_change,prct_change,...
             prct,prct,prct,prct};


nvar=length(brsel_all);   

brseries_gr=zeros(N_economy, NT_sim, nvar);
for s=1:N_economy
    for v=1:nvar
        thisv=brsel_all(v); 
        this_series = simseries_mean_econ{s}(1:NT_sim,thisv);
        this_series = transformvec{v}(this_series);
        brseries_gr(s,:,v) = this_series;
    end
end

% add data
datatab=readtable(datapath);
datavals=datatab{3:NT_sim+2,{'sbankds','sbanklev','cbanklev'}};
datavals=datavals*100;
datavals=[nan(NT_sim,3), datavals];

brseries_all = zeros(N_economy+1, NT_sim, nvar);
brseries_all(1,:,:) = datavals;
brseries_all(2:end,:,:) = brseries_gr;

titles4={'BDS Output', 'Consumption','Liquidity', ...
         'S Debt Share','S Leverage', 'C Leverage'}; 

     
if usejava('desktop')
    makeImpRes(brseries_all,tvec,titles4,colors,[2,3],label,printfiles{1});   
end

