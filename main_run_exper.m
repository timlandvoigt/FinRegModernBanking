if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear;
end
close all;

% ===========================================
% program control
% ===========================================

% path to file with experiment definition
if ~exist('exper_path','var')
    exper_path='env_base.mat';
end
if ~exist('maxit','var')
    maxit=10;
end
if ~exist('tol_avg','var')
    % mean convergence
    tol_avg=1e-5;
end
if ~exist('localpool','var')
    % local cluster profile for parpool?
    localpool=true;
end
if ~exist('no_par_processes','var')
    % local cluster profile for parpool?
    no_par_processes=4;
end


open_parpool; 

% print mode
printmode=1;

% ===========================================
% load environment structure
% ===========================================

load(exper_path);

% ===========================================
% outer convergence loop
% ===========================================

% policy iteration
mobj=mobj.polIter(maxit,1,printmode,tol_avg);


% ===========================================
% results
% ===========================================


% make date string
curr_date=clock;
date_str='';
for i=1:5
   date_str=[date_str,'_',num2str(curr_date(i))]; 
end

if ~exist('outname','var')
    outname=['res',date_str,'.mat'];
else
    outname=[outname,date_str,'.mat'];    
end

save(outname,'mobj','stv');              




