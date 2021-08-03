if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear;
end

% ===========================================
% program control
% ===========================================

% name of file with experiment definitions
if ~exist('experdef_file','var')
        experdef_file='experdef_20201117.m';
end
% name of experiment
if ~exist('expername','var')
             expername='base'; 
end

% possible values 'no_guess', 'guess'
if ~exist('guess_mode','var')
          guess_mode='no_guess';
end
% path to file with initial guess; not used if guess_mode='no_guess'
if ~exist('guess_path','var')
     guess_path = '';
end

% path and file name for output file
outfname=['env_',expername,'.mat'];

% approximation mode
approxMode='linear';

% only compute steady state (without writing initial model object)
ststonly=0;

run(experdef_file)
disp(experdef_file);
expdef=allexpers.(expername);
params=expdef.params;

if ~params.C_bank_only
    modelclassname='SBModel_CD'; 
else
    modelclassname='CModel'; 
end
gvec = expdef.GuessVec;   
if params.calib_eps
    eps=0.3;
    gvec=[gvec, -log((1-eps)/(1+eps))];
end
ststfun=str2func([modelclassname,'.compStSt']);

instfun=str2func(modelclassname);
guessfun=str2func([modelclassname,'.assignGuess']);

% compute steady state values
% do this once for each state of sigma_omega
options=optimset('Display','iter','TolX',1e-10,'TolFun',1e-10, ...
    'MaxIter',3e3, 'MaxFunEvals',3e3);


for rstate =   1 
    
    expdef.params.varrhoS_ss = params.varrho_val(rstate) ; % steady state value
    
    disp('----------------------------------------------------')
    disp('----------------------------------------------------')
    disp([num2str(expdef.params.varrhoS_ss,2) ' % of hh run'])
    disp('----------------------------------------------------')
    disp('----------------------------------------------------')
    
    fh_compStSt=@(x)ststfun(x,expdef.params,0);
    [sol_here,~,exfl]=fsolve(fh_compStSt,gvec,options);
    if exfl<1
        disp('!! Problem computing steady state');
    end
    [~,stv]=ststfun(sol_here,expdef.params,2);
    
end



if ststonly
    return;
end


% ===================================================
% Parameters of stochastic model
% ===================================================

% for GDP component
N_y=3;
sig2_y=params.sigY^2;
mu_Y=params.muY ;
rho_Y=params.delta_Y ;
% sig2_y=log(sig2_Y/mu_Y +1);
mu_y=log(mu_Y)-0.5*sig2_y;
Skew_y=0;
[Yprob,y] = model.DSGEModel.rouwen(rho_Y,mu_y,sqrt(sig2_y),Skew_y,N_y);
Yprob=Yprob';
Y=exp(y);

% for Z shock to capital 
N_z=3;
[Zprob,z] = model.DSGEModel.rouwen(0,0,params.sigZ,Skew_y,N_z);
Z = exp(z);
% time-varying cap req
vartheta = params.vartheta;
vartheta_vec=[vartheta(1); sum(vartheta)/2; vartheta(2)];

% for Rho shock to run risk
N_rho=2;
runprob = params.varrhoprob;
varrho = params.varrho_val ;

% Markov transition for Y and Z
comb=grid.StateSpaceGrid.makeCombinations([N_y,N_z]);
comstates=[ Y(comb(:,1)), Z(comb(:,2)), vartheta_vec(comb(:,2)) ]; % cap req perfectly correl. with Z
N_comb=N_y*N_z;
comtrans=kron(Yprob ,Zprob);
trans_sum = repmat(sum(comtrans,2),1,size(comtrans,2));
comtrans = comtrans./trans_sum ;

% bank runs only in bad payoff states for Z
bust_ind=comstates(:,2)<1;
N_bust=sum(bust_ind);
bust_comb=comstates(bust_ind,:);
boom_comb=comstates(~bust_ind,:);
N_boom=N_comb-N_bust;
rhocomb=grid.StateSpaceGrid.makeCombinations([N_bust,N_rho]);
rhostates=[ bust_comb(rhocomb(:,1),:), varrho(rhocomb(:,2))];
boomstates=[ boom_comb, ones(N_boom,1)*varrho(1) ];
mpts_perm=[ rhostates; boomstates];
mpts_perm=[ mpts_perm(:,1:2), mpts_perm(:,4), mpts_perm(:,3) ];

% transition matrix
rhotrans= [ kron(comtrans(bust_ind,bust_ind),runprob), kron(comtrans(bust_ind,~bust_ind),ones(N_rho,1)) ];
boomtrans=[ kron(comtrans(~bust_ind,bust_ind), runprob(1,:)), comtrans(~bust_ind, ~bust_ind)];
mtrans=[rhotrans; boomtrans];


% correlation structure between Z and Y
exnpt=size(mpts_perm,1);
mpts_perm(:,2)=params.scaleZ*mpts_perm(:,2).*mpts_perm(:,1);
mpts_all=[mpts_perm, (1:exnpt)'];
     

% simulate exog process and produce histogram
Nsim=10000;
simst=model.DSGEModel.hitm_s(mtrans',rand(Nsim,1));
simfrac=histc(simst,1:exnpt)/Nsim;
disp('-----------------------------------------');
disp('Unconditional prob of each state: ');
disp(num2str([(1:exnpt)',mpts_perm,simfrac]));

% variable lists
% exogenous state variable names
exogenv=struct;
exogenv.exnames={'Y','Z','varrho','vartheta'};
exogenv.exnpt=exnpt;
exogenv.pts_perm=mpts_perm;
exogenv.pts_all=mpts_all;
exogenv.mtrans=mtrans;

          
% ===========================================
% create model object
% ===========================================

% endogenous state variables
endogenv=struct;
% assign values for initial guess
[solguessvec,Vguessvec,Vnames]=guessfun(stv);

% save name lists
endogenv.solnames=fieldnames(stv.Sol);
endogenv.solbase=solguessvec;
endogenv.addnames=fieldnames(stv.Add);
endogenv.Vnames=Vnames;

if params.C_bank_only 
    endogenv.ennames={'K','bC'};
    unigrids={1:exogenv.exnpt,expdef.Kpts,expdef.bCpts};    
    
else
    endogenv.ennames={'K','bS','bC','KSsh'};
    unigrids={1:exogenv.exnpt,expdef.Kpts,expdef.bSpts,expdef.bCpts,expdef.KSshpts};
end
basegrid=grid.TensorGrid(unigrids);
mgrid=basegrid;

if isequal(guess_mode,'no_guess')
    solmatguess=repmat(solguessvec',[mgrid.Npt,1]);
    Vmatguess=repmat(Vguessvec',[mgrid.Npt,1]);
    
else
    guessobj=load(guess_path,'mobj');
    solmatguess=guessobj.mobj.evaluatePol(mgrid.Pointmat)';
    Vmatguess=guessobj.mobj.evaluateVal(mgrid.Pointmat)';
end
            
             
% build approximating function
% create guess for solution variable functions
Pf=grid.LinearInterpFunction(mgrid,solmatguess,true);
% create guess for next-period functions
Vf=grid.LinearInterpFunction(mgrid,Vmatguess,true);
% and for state transition function
Tf=[];
            
mobj=instfun(expdef.params,endogenv,exogenv,Vf,Pf,Tf);
disp('-------------------------------------');
disp('Bounds of state space:');
disp(num2str(mobj.Vfct.SSGrid.StateBounds));

% ===========================================
% save initial object
% ===========================================

save(outfname,'stv','mobj');              






