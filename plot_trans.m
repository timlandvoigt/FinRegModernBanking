close all;
clear;

respath='./';
outpath='./Results/';
resfile='res_20201117_base';
outfile=['GTR_',resfile];

load([respath,resfile,'.mat']);
load([respath,'GTR_',resfile,'.mat']);


outpath=[outpath,outfile,'_'];
tvec=0:NT_sim-1;



%% Make other Graphs

% file names for graphs (set to empty for no printing)
printfiles={[outpath,'IRF1'],[outpath,'IRF2'],[outpath,'IRF3'],[outpath,'IRF4'],[outpath,'IRF5']};        

% which variables
brsel1=[indexmap.get('Z'),indexmap.get('I'),indexmap.get('GDP'),...
        indexmap.get('C'),indexmap.get('H'),indexmap.get('DWL')];  
brsel2=[indexmap.get('KSsh'),indexmap.get('Slev'),indexmap.get('Clev'),...
        indexmap.get('FS'),indexmap.get('FC'),indexmap.get('rateS')];  
brsel_all=[brsel1];
nvar=length(brsel_all);   

relative_plot=[ones(1,5),2];

% How many shocks to plot (one less for relative)
N_shock_plot = N_shock - 1;

% transforms
plain = @(x)(x);
prct = @(x)(100*x);
prct_change = @(x)((x/x(1) -1)*100);
transformvec={plain,plain,plain,plain,plain,plain};


nvar=length(brsel_all);   

brseries_gr=zeros(N_shock_plot, NT_sim, nvar);
for s=1:N_shock_plot
    for v=1:nvar
        thisv=brsel_all(v); 
        if relative_plot(v)==1
            this_series = simseries_diff_mean{s+1}(1:NT_sim,thisv)./simseries_mean{1}(1:NT_sim,thisv) * 100;
        elseif relative_plot(v)==2
            this_series = simseries_diff_mean{s+1}(1:NT_sim,thisv) * 100;
        else
            this_series = simseries_mean{s}(1:NT_sim,thisv);
        end
        this_series = transformvec{v}(this_series);
        brseries_gr(s,:,v) = this_series;
    end
end


colors={'k-o','r-o'};
if N_shock_plot==3
   colors = ['b-o',colors]; 
end

titles1={'Z', 'BDS Investment', 'BDS Output', 'Consumption',  'Liquidity', 'DWL/GDP'}; 
titles2={'S Capital Share', 'S Leverage', 'C Leverage', 'S Default Rate', 'C Default Rate', 'S Deposit Rate'}; 

if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:6),tvec,titles1,colors,[2,3],[],printfiles{1});   
    %makeImpRes(brseries_gr(:,:,7:12),tvec,titles2,colors,[2,3],[],printfiles{2});  
end

