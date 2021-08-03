% set params structure
params=struct;

% --- Model select 
params.use_kappa_fair=false;
params.C_bank_only=false;
params.CD_U    = true;
params.use_theta_var=false;
params.Splanner = false;
params.norun    = false ;
params.highshocks = false;
params.calib_eps = false;
params.betaScoef = -0.0019 ; % coefficient from data regression

% --- Benchmark welfare
params.benchVH = -141.3613;

params.muY      = 1      ; % normalization
if params.highshocks
params.sigY     =  0.0087*5  ; % backed out std from GDP 1999-present
else
params.sigY     =  0.0087  ; % backed out std from GDP 1999-present
end
params.delta_Y  =  0.5988  ; % AC

params.scaleZ   =  0.196  ; % Z*K/GDP ... share of fin value added
params.muZ      = params.muY*params.scaleZ  ;
if params.highshocks
params.sigZ     =0.0187*5 ; % backed out std from bank dep sales growth 1999-present
else
params.sigZ     =0.0174; % backed out std from bank dep sales growth 1999-present
end

params.Zbar_disc = .26 ; 
params.muZbar=params.muZ*params.Zbar_disc;
if params.norun
    params.varrho_val = [0;0];
else
   params.varrho_val = [0;1/3];
end
params.varrhoprob = [.9625 1-.9625; %  0.0476% unconditional run prop ... 
                     3/4 1/4]; % average run lasts 4 months
% params.varrhoprob = [0.95 0.05; % 5% disaster prop ... 
%                      0.33  0.67]; % probability of run65
params.deltaKbar = 0.147;
params.thetaS    = 0.99 ;
params.mu_rho   = [1,1];

params.gamma  =  2   ;       % bob hall parameter
params.gammaH =  1.6 ;% 1.6   ;     % bob hall parameter
params.eta    = 2/3  ;       % labor share
params.theta  = 0.1  ;       % cap req
params.vartheta=repmat(params.theta,1,2);
params.kappaC  = 0.00142 ;  % assessement fee
params.kappaS  = 0.0 ;  % tax on S bank debt
params.FCbar  = 0.001      ; % target default rate
params.deltaK = 0.025; % 10% annual depreciation

% -- calibration targets
params.beta   = 0.993  ;%0.9919 ; % 0.989    ; % get risk free rate:  average quarterly Tbill 99-15
params.psi    = 0.0072; %.008   ; % liquidity prem
params.nu     =  0; % weight on  
params.fullnu =  0; % also have run and bailout prob in liquidity factor
params.piB    = .85       ; %   63  0.9809
% params.piB    = 0.8  ; %   63  0.9809
params.xiC    = 0.352   ;    % recovery rate depository ~60/70%
params.xiS    = 0.2051   ;% 0.38  ;    % recovery rate Moody financial firms (~ 37%)
%params.xiS    = 0.05  ;    % recovery rate Moody financial firms (~ 37%)
params.deltaS  = 0.39     ;   % default rate S 
% params.deltaS  = 0   ;   % default rate S 
params.deltaC  = 0.204  ; % default rate C 
params.alpha   = 0.33   ; % share shadow bank activity
params.epsilon = 0.2   ; 
params.Nbar = 1;

% -- dynamic targets.. guess 
params.phi     = 0.3  ; % 2nd moment - investment vol
params.phiK    = 0.011  ; % 025 
params.sig_rho = [ 0.121    ....
                   0.254 ];
              
% Steady state 
GuessVec =[0.0015   -0.1227    0.2795    1.1771    0.4341   -1.6762    1.6762   -1.2039   -1.5327];
GuessVec_SP=[.9 .9 1 1 0.5];
             
% --- other definition
params.mu_rhoC  = params.mu_rho(1);
params.mu_rhoS  = params.mu_rho(2);
params.sig_rhoC = params.sig_rho(1);
params.sig_rhoS = params.sig_rho(2);

basegrid.Kpts = linspace(2.92,3.5   ,6); % gammaH  =1.5
basegrid.bSpts = linspace(0.1,.92,8);
basegrid.bCpts = linspace(0.87,.93,6);
basegrid.KSshpts = linspace(0.2,0.42 ,7);    



% base experiment definition
xpbase=struct('params',params,...
              'Kpts',basegrid.Kpts,...
              'bSpts',basegrid.bSpts,...
              'bCpts',basegrid.bCpts,...
              'KSshpts',basegrid.KSshpts,...
              'GuessVec',GuessVec);
          
% baseline theta variations
xptheta13=xpbase;
xptheta13.params.theta=0.13;
xptheta13.bCpts = linspace(0.85,.89,6);

xptheta14=xpbase;
xptheta14.params.theta=0.14;
xptheta14.bCpts = linspace(0.83,.89,6);

xptheta15=xpbase;
xptheta15.params.theta=0.15;
xptheta15.bCpts = linspace(0.82,.88,6);

xptheta16=xptheta15;
xptheta16.params.theta=0.16;
xptheta16.bCpts = linspace(0.82,.87,6);

xptheta17=xptheta15;
xptheta17.params.theta=0.17;
xptheta17.bCpts = linspace(0.81,.85,6);

xptheta18=xptheta15;
xptheta18.params.theta=0.18;
xptheta18.bCpts = linspace(0.8,.84,6);

xptheta20=xpbase;
xptheta20.params.theta=0.2;
xptheta20.bCpts = linspace(0.77 ,.83,6);

xptheta25=xpbase;
xptheta25.params.theta=0.25;
xptheta25.Kpts = linspace(3.00,3.5,6);    
xptheta25.bSpts =  linspace(0.1,0.96,8);
xptheta25.bCpts = linspace(0.74 ,.82,6);

xptheta30=xptheta25;
xptheta30.params.theta=0.3;
xptheta30.bCpts = linspace(0.69 ,.73,6);

% gammaH=0
xpbasegH0=xpbase;
xpbasegH0.params.gammaH=0;
xpbasegH0.params.psi=0.004;
xpbasegH0.Kpts = linspace(2.9,3.3,6);
xpbasegH0.bSpts = linspace(0.04,0.9,6);

xpbasegH0.GuessVec=[0.0002   -0.1128    0.1885    0.6680    0.0841   -1.2191    1.2191   -1.0678   -1.2940];


xptheta13gH0=xpbasegH0;
xptheta13gH0.params.theta=0.13;
xptheta13gH0.bCpts= linspace(0.85,.89,6);

xptheta14gH0=xpbasegH0;
xptheta14gH0.params.theta=0.14;
xptheta14gH0.bCpts= linspace(0.84,.88,6);


xptheta15gH0=xpbasegH0;
xptheta15gH0.params.theta=0.15;
xptheta15gH0.bCpts= linspace(0.83,.88,6);

xptheta16gH0=xpbasegH0;
xptheta16gH0.params.theta=0.16;
xptheta16gH0.bCpts= linspace(0.82,.87,6);


% C-bank only
xpbaseC=xpbase;
xpbaseC.params.C_bank_only = true;
xpbaseC.params.gammaH=0;
xpbaseC.params.psi=0.00225;
xpbaseC.Kpts= linspace(2.7,3.3,7);
xpbaseC.GuessVec = [log(1), log(0.98), log(3*params.theta), 0.5];
xpbaseC.bCpts= linspace(0.87,.92,6);
xpbaseC.Kpts= linspace(3.00,3.5,6);

xptheta13C=xpbaseC;
xptheta13C.params.theta=0.13;
xptheta13C.bCpts= linspace(0.85,.92,6);

xptheta14C=xpbaseC;
xptheta14C.params.theta=0.14;
xptheta14C.bCpts= linspace(0.84,.91,6);

xptheta15C=xpbaseC;
xptheta15C.params.theta=0.15;
xptheta15C.bCpts= linspace(0.83,.9,6);

xptheta16C=xpbaseC;
xptheta16C.params.theta=0.16;
xptheta16C.bCpts= linspace(0.83,.89,6);


% Precrisis
xpprecrisis=xpbase;
xpprecrisis.params.piB = 0.87;
xpprecrisis.params.deltaS = 0.45;
xpprecrisis.params.theta = 0.08;
xpprecrisis.params.varrho_val = [0;0];
xpprecrisis.bCpts = linspace(0.89,.95,6);
xpprecrisis.bSpts = linspace(0.2,.95,6);
xpprecrisis.GuessVec = [0.0003   -0.1157    0.1865    0.5679    0.1275   -1.1544    1.1544   -1.0844   -1.2378];

xppostcrisis=xpbase;
xppostcrisis.params.piB = 0.7;
xppostcrisis.params.theta = 0.11;
xppostcrisis.bCpts = linspace(0.85,.9,6);
xppostcrisis.bSpts = linspace(0.2,.85,6);

xppostcrisis08=xppostcrisis;
xppostcrisis08.params.theta = 0.08;
xppostcrisis08.params.piB = 0.7;
xppostcrisis08.bCpts = linspace(0.89,.95,6);
xppostcrisis08.bSpts = linspace(0.2,.85,6);


xppostcrisis085=xppostcrisis;
xppostcrisis085.params.theta = 0.085;
xppostcrisis085.params.piB = 0.84;
xppostcrisis085.bCpts = linspace(0.89,.95,6);
xppostcrisis085.bSpts = linspace(0.2,.93,6);

xppostcrisis09=xppostcrisis;
xppostcrisis09.params.theta = 0.09;
xppostcrisis09.bCpts = linspace(0.88,.94,6);
xppostcrisis09.params.piB = 0.81;
xppostcrisis09.bSpts = linspace(0.2,.9,6);

xppostcrisis095=xppostcrisis;
xppostcrisis095.params.theta = 0.095;
xppostcrisis095.bCpts = linspace(0.87,.93,6);
xppostcrisis095.params.piB = 0.77;
xppostcrisis095.bSpts = linspace(0.2,.89,6);

xppostcrisis10=xppostcrisis;
xppostcrisis10.params.theta = 0.1;
xppostcrisis10.bCpts = linspace(0.87,.92,6);
xppostcrisis10.params.piB = 0.74;
xppostcrisis10.bSpts = linspace(0.2,.89,6);

xppostcrisis105=xppostcrisis;
xppostcrisis105.params.theta = 0.105;
xppostcrisis105.bCpts = linspace(0.87,.92,6);
xppostcrisis105.params.piB = 0.71;
xppostcrisis105.bSpts = linspace(0.2,.89,6);

xppostcrisisc1=xppostcrisis;
xppostcrisisc1.params.theta = 0.08;
xppostcrisisc1.bCpts = linspace(0.89,.95,6);
xppostcrisisc1.params.piB = 0.84;

xppostcrisisc2=xppostcrisisc1;
xppostcrisisc2.params.piB = 0.81;

xppostcrisisc3=xppostcrisisc1;
xppostcrisisc3.params.piB = 0.77;

xppostcrisisc4=xppostcrisisc1;
xppostcrisisc4.params.piB = 0.74;

xppostcrisisc5=xppostcrisisc1;
xppostcrisisc5.params.piB = 0.71;
%xppostcrisis.GuessVec = [0.0003   -0.1307    0.1867    0.5664    0.1165   -1.1551    1.1551   -1.0859   -1.2372];



% simple: equal xi, delta and rho + piB = 0 and no runs
% both banks like C bank
xpbasesimpleCgH0=xpbase;
xpbasesimpleCgH0.params.xiS = xpbase.params.xiC;
xpbasesimpleCgH0.params.deltaS = xpbase.params.deltaC;
xpbasesimpleCgH0.params.sig_rhoS=xpbase.params.sig_rhoC;
xpbasesimpleCgH0.params.piB = 0;
xpbasesimpleCgH0.params.gammaH=0;
xpbasesimpleCgH0.params.psi=0.004;
xpbasesimpleCgH0.params.varrho_val = [0;0];
xpbasesimpleCgH0.bSpts= linspace(0.65,.83,6);
xpbasesimpleCgH0.KSshpts= linspace(0.35,.45,7);

xptheta13simpleCgH0=xpbasesimpleCgH0;
xptheta13simpleCgH0.params.theta = 0.13;
xptheta13simpleCgH0.bCpts= linspace(0.85,.89,6);

xptheta14simpleCgH0=xpbasesimpleCgH0;
xptheta14simpleCgH0.params.theta = 0.14;
xptheta14simpleCgH0.bCpts= linspace(0.84,.88,6);


xptheta15simpleCgH0=xpbasesimpleCgH0;
xptheta15simpleCgH0.params.theta = 0.15;
xptheta15simpleCgH0.bCpts= linspace(0.83,.88,6);

xptheta16simpleCgH0=xpbasesimpleCgH0;
xptheta16simpleCgH0.params.theta = 0.16;
xptheta16simpleCgH0.bCpts= linspace(0.82,.86,6);



% nu = 1
xpbasenu=xpbase;
xpbasenu.params.nu = 1;

% nu = 10
xpnu10=xpbase;
xpnu10.params.nu = 10;

% fullnu
xpfullnu2=xpbase;
xpfullnu2.params.nu=2;
xpfullnu2.params.fullnu=1;


xppiBup=xpbase;
xppiBup.params.piB=.865;
xppiBup.bSpts = linspace(0.3,0.93,7);
xppiBup.KSshpts = linspace(0.25,0.42,7);


xppiBdn=xpbase;
xppiBdn.params.piB=0;
xppiBdn.bSpts = linspace(0.3,0.8,7);


% structure with all experiments
allexpers=struct('base',xpbase,...
                 'baseC',xpbaseC,...
                 'theta14C',xptheta14C,...
                 'theta15C',xptheta15C,...
                 'theta16C',xptheta16C,...
                 'theta13',xptheta13,...
                 'theta14',xptheta14,...
                 'theta15',xptheta15,...
                 'theta16',xptheta16,...
                 'theta17',xptheta17,...
                 'theta20',xptheta20,...
                 'theta30',xptheta30,...
                 'basegH0',xpbasegH0,...
                 'theta14gH0',xptheta14gH0,...
                 'theta15gH0',xptheta15gH0,...
                 'theta16gH0',xptheta16gH0,...
                 'piBup',xppiBup,...
                 'piBdn',xppiBdn,...
                 'basenu',xpbasenu,...
                 'nu10',xpnu10,...
                 'fullnu2',xpfullnu2,...
                 'basesimpleCgH0',xpbasesimpleCgH0,...
                 'theta14simpleCgH0',xptheta14simpleCgH0,...
                 'theta15simpleCgH0',xptheta15simpleCgH0,...
                 'theta16simpleCgH0',xptheta16simpleCgH0,...
                 'precrisis',xpprecrisis,...
                 'postcrisis',xppostcrisis,...
                 'postcrisis08',xppostcrisis,...
                 'postcrisis085',xppostcrisis085,...
                 'postcrisis09',xppostcrisis09,...
                 'postcrisis095',xppostcrisis095,...
                 'postcrisis10',xppostcrisis10,...
                 'postcrisis105',xppostcrisis105,...
                 'postcrisisc1',xppostcrisisc1,...
                 'postcrisisc2',xppostcrisisc2,...
                 'postcrisisc3',xppostcrisisc3,...
                 'postcrisisc4',xppostcrisisc4,...
                 'postcrisisc5',xppostcrisisc5);

% % parameter sensitivity
sscnames={'sig_rhoS'};
sscvar=0.05;
updn={'up','dn'};
updnsgn={1,-1};
sscexpers=cell(1,numel(sscnames)*2);
ssccnt=1;
for i=1:numel(sscnames)
    thisparname=sscnames{i};
    thisexpname=thisparname;
    thisexpname(regexp(thisexpname,'[_]'))=[];
    for d=1:2
       thisexp=xpbase;
       thisexp.params.(thisparname)=xpbase.params.(thisparname)*(1+sscvar*updnsgn{d});
       thisname=['ssc',thisexpname,updn{d}];
       allexpers.(thisname)=thisexp;
       sscexpers{ssccnt}=thisname;
       ssccnt=ssccnt+1;
    end
end


% custom grid for sensitivity checks            
% experselect={'theta13','theta14','theta15','theta16','theta17','theta20','theta30',...
%              'postcrisisc1','postcrisisc2','postcrisisc3','postcrisisc4','postcrisisc5'};
% 
% experselect=[{'piBup','piBdn'}, sscexpers];
      
% experselect={'eps1','eps2','eps3','eps4',...
%               'gamH1','gamH2','gamH3','gamH4'};

%experselect={'baseC','theta13C','theta14C','theta15C','theta14gH0','theta14simpleCgH0'};
%experselect=sscexpers;
%experselect=[{'minn16'}, sscexpers];
%experselect={'basexidelta','theta15xidelta','theta20xidelta','theta25xidelta','theta30xidelta'};
%experselect={'basesimpleCgH0','theta13simpleCgH0','theta15simpleCgH0','theta16simpleCgH0','theta20simpleCgH0','theta25simpleCgH0'};
%experselect={'basegH0','theta13gH0','theta15gH0','theta20gH0','theta25gH0'};
   %experselect={'theta13piB','theta20simple','theta25simple','theta30simple'};

         
expernames=fieldnames(allexpers);
experselect = expernames;

% Write list of defined experiments to file
fid = fopen([mfilename,'.txt'],'w');
for i=1:length(experselect)
   fprintf(fid,experselect{i});
   fprintf(fid,'\n');
end
fclose(fid);             
          
          



