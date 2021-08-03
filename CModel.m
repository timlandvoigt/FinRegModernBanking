
classdef CModel < model.DSGEModel
%     
    properties (SetAccess=protected)
        NSTEN % number of endogenous state variables: Vfct.Ndim-1 
        NSTEX % number of exogenous state variables
        NSOL % number of solution vars
        NV % number of forecasting variables
        NADD % number of additional endogenous variables
        Sol_names % NSOLx1 cell array with names of solution variables
                   % must be in order of functions of Pfct 
        V_names % NVx1 cell array with names of forecasting variables
                   % must be in order of functions of Vfct      
        Sol_baseguess           
        En_names % NSTENx1 cell array with names of endog state variables           
        Ex_names % NSTEXx1 cell array with names of exog state variables           
        Add_names % NADDx1 cell array with names of additional vars
        Params % structure array with parameters
        Exogenv % structure array with fields mtrans, pts_perm, pts_all
                % specifying exog. Markov processes
        Vfct % ApproxFunction object for iteration-relevant functions
        Pfct % ApproxFunction object for solution jump variables
        Tfct % ApproxFunction object for transition of state variable(s) (optional)
    end    
    
    
    methods
        % constructor
        function obj=CModel(params,endogenv,exogenv,vfct,pfct,tfct)
            % call superclass constructor
            obj=obj@model.DSGEModel(params,endogenv,exogenv,vfct,pfct,tfct);
        end
        
        function [nextst,outstr]=calcStateTransition(obj,point,solvec,mode,varargin)

            % unpack params
            params = obj.Params;
            phi       = params.phi;
            phiK      = params.phiK;
            nu        = params.nu;
            alpha     = params.alpha;
            eta       = params.eta  ; 
            deltaK    = params.deltaK;
            deltaKbar = params.deltaKbar;
            if isfield(params,'kappaC')
                kappaC     = params.kappaC ;
                kappaS     = params.kappaS ;
            else
                kappaC=params.kappa;
                kappaS=0;
            end
            epsilon   = params.epsilon;
            xiC       = params.xiC ;
            xiS       = params.xiS ;
            Zbar_disc = params.Zbar_disc;
            mu_rhoC   = params.mu_rhoC ;
            sig_rhoC  = params.sig_rhoC ;
            mu_rhoS   = params.mu_rhoS ;
            sig_rhoS  = params.sig_rhoS ;
            deltaS    = params.deltaS;
            deltaC    = params.deltaC;
            chi1C     = sig_rhoC^2/mu_rhoC;
            chi0C     = mu_rhoC/chi1C;
            chi1S     = sig_rhoS^2/mu_rhoS;
            chi0S     = mu_rhoS/chi1S ;
            Nbar = params.Nbar;
             % extract state variables
            exst=point(1);
            K=point(2);
            bC=point(3);
                        
            Y=obj.Exogenv.pts_perm(exst,1);
            Z=obj.Exogenv.pts_perm(exst,2);
            Zbar=Zbar_disc*Z;
            varrho=obj.Exogenv.pts_perm(exst,3);
            if params.use_theta_var
                vartheta=obj.Exogenv.pts_perm(exst,4);
            else
                vartheta=params.theta;
            end
            
            if isempty(varargin)
                thisvals=obj.evaluateVal(point);
            else
                thisvals=varargin{1};
            end        
            
            % extract solution variables
            p = exp(solvec(1));
            qC = exp(solvec(2));
            BCpol = exp(solvec(3));
                        
            % small value functions for default cutoffs
            vC=phiK/2*(thisvals(5)^2-1);
             
            % this period
            KC = K;
            NC = Nbar;
            nC = NC/KC;
            
            %nH = nS*(Z/Zbar)^(1/(eta-1));
            
%             PiS = (1-eta)*Z*nS^eta  + p - deltaK + ((p-1)^2)/(2*phi) ; 
%              LS = bS/PiS ;
%             PiH = (1-eta)*Zbar*nH^eta  + p *(1-deltaKbar) ; 
%              xt = PiH./PiS;
%              lS = varrho.*LS./xt ;
%              NS = nS*KS*(1-lS) ;
%              NH = nH*lS*KS ;
%          DFcutS = ((1-varrho).*LS - (1-lS)*(vS/PiS + deltaS))./(1-lS) ;   % REV: new default cutoff
%              FS = gamcdf(DFcutS,chi0S,chi1S);
%            FSrho_minusS = mu_rhoS*gamcdf(DFcutS,chi0S+1,chi1S);
%            FSrho_plusS  = mu_rhoS*(1-gamcdf(DFcutS,chi0S+1,chi1S));
%            
             
            PiC = (1-eta)*Z*nC^eta  + p - deltaK + ((p-1)^2)/(2*phi) ; 
             BC = bC*KC;
             LC = bC/PiC;
         DFcutC = LC-deltaC - vC/PiC; 
             FC = gamcdf(DFcutC ,chi0C,chi1C);
   FCrho_minusC = mu_rhoC*gamcdf(DFcutC ,chi0C+1,chi1C);
    FCrho_plusC = mu_rhoC*(1-gamcdf(DFcutC ,chi0C+1,chi1C));
            
            % investment 
            IC  = ((p-1)/phi  + deltaK)*KC ;
            i_c = IC/KC ;

            % deadweight losses 
            DWL_C = xiC*FCrho_minusC* (PiC-(1-deltaK)*p)*KC ;
           % capital transition
            Kpol =   IC  + (1-deltaK)*(1-xiC*FCrho_minusC)*KC;
           KCpol = Kpol;
           kC = KCpol/KC;
           
           % consumption
           phiI_C = (phi/2)*((i_c - deltaK)^2)*KC +        phiK/2*(kC-1)^2 *KC        ;
           Y_C = Z*(KC^(1-eta))*NC^eta ;
           C    = Y    +  Y_C    +  ....
               - IC     - phiI_C ...
               - DWL_C ;
           % liquidity
           
%            if epsilon==0
%                H = (BS.^alpha_eff).*(BC^(1-alpha_eff));
%            else
%                H = (alpha_eff*BS^epsilon + (1-alpha_eff)*BC^epsilon)^(1/epsilon) ;
%            end
           H = BC;
           
           
           % dividends
           DC = FCrho_plusC*KC*PiC - (1-FC)*BC + (1-FC)*((qC-kappaC)*BCpol - p*KCpol - phiK/2*(kC-1)^2 *KC );
           
           
           exnpt = obj.Exogenv.exnpt;

            if mode>0
    % simulation, mode contains number of next period's state
                cind=obj.Exogenv.pts_all(mode,end);
                Knext=Kpol;
                bCnext=BCpol/KCpol;
            else
                % solution mode, compute next period's state variables for all
                % possible Markov states
                cind=obj.Exogenv.pts_all(:,end);
                Knext=Kpol*ones(exnpt,1);
                bCnext=BCpol/KCpol*ones(exnpt,1);
            end
            
            nextst=[cind,Knext,bCnext];            
            
            addvars=struct('Knext',Kpol,...
                'BCnext',BCpol,...
                'kC',kC,...
                'H',H,...
                'C',C,...
                'DC',DC,...
                'Z',Z,...
                'Y',Y,....
                'NC',NC);
            
            outstr=struct;
            outstr.addvars=addvars;
            outstr.exstvec=[Y;Z;varrho;vartheta];

        end
        

        
        function [fx,V]=calcEquations(obj,exst,nextst,solvec,instr,mode)
            
            % allocate result
            fx=zeros(obj.NSOL,1);
            
            % unpack params
            params=obj.Params;
            alpha     = params.alpha;
            psi       = params.psi   ; % weight on liquidity
            beta      = params.beta ;
            eta       = params.eta;
            Nbar      = params.Nbar;
            nu        = params.nu;
            phi       = params.phi;
            phiK      = params.phiK;
            deltaK    = params.deltaK;
            deltaKbar = params.deltaKbar;
            Zbar_disc = params.Zbar_disc;
            if isfield(params,'kappaC')
                kappaC     = params.kappaC ;
                kappaS     = params.kappaS ;
            else
                kappaC=params.kappa;
                kappaS=0;
            end
            theta     = params.theta ;
            thetaS    = params.thetaS;
            gamma     = params.gamma ;
            gammaH    = params.gammaH ;
            xiC       = params.xiC ;
            xiS       = params.xiS ;
            mu_rhoC   = params.mu_rhoC ;
            sig_rhoC  = params.sig_rhoC ;
            mu_rhoS   = params.mu_rhoS ;
            sig_rhoS  = params.sig_rhoS ;
            deltaS    = params.deltaS;
            deltaC    = params.deltaC;
            piB       = params.piB;
            if isfield(params,'epsilon')
                epsilon=params.epsilon;
            else
                epsilon=0;
            end
            chi1C     = sig_rhoC^2/mu_rhoC;
            chi0C     = mu_rhoC/chi1C;
            chi1S     = sig_rhoS^2/mu_rhoS;
            chi0S     = mu_rhoS/chi1S;
            % extract endogeous variables
            p    = exp(solvec(1));
            qC   = exp(solvec(2));
            lamC = solvec(4);
            
            lamCplus=max(0,lamC)^3;
            lamCminus=max(0,-lamC)^3;
            
            
            % extract some state-dependent values
            envec = instr.addvars;
            Knext=envec.Knext;
            BCnext = envec.BCnext;
            C = envec.C;
            H = envec.H;
            NC = envec.NC;
            KCnext = Knext;
            kC = envec.kC;
            if params.use_theta_var
                theta=instr.exstvec(4);
            end
            
            % compute expectation terms
            prnext = obj.Exogenv.mtrans(exst,:);
            Ynext  = obj.Exogenv.pts_perm(:,1) ;
            Znext  = obj.Exogenv.pts_perm(:,2) ;
            
            %Zbarnext=obj.Exogenv.pts_perm(:,3);
%             Zbarnext   = Zbar_disc*Znext;
%             varrhonext = obj.Exogenv.pts_perm(:,3);
            
            % projection evaluation
            Pol_next   = obj.evaluateVal(nextst)';
            pnext      = Pol_next(:,1);
            kCnext     = Pol_next(:,5);
            vCnext = phiK/2*(kCnext.^2 - 1);
                        
            NCnext = Nbar*ones(size(Znext));
            nCnext = NCnext ./ KCnext;
            
            % next period
            ICnext  = ((pnext-1)/phi  + deltaK).*KCnext ;
            i_cnext = ICnext/KCnext ;
            
%             DFcutS     = ((1-varrhonext).*LSnext-(1-lSnext).*(vSnext./PiSnext+deltaS))./(1-lSnext) ;  % REV: new default threshold
%             FSnext     = gamcdf(DFcutS,chi0S,chi1S);
%             fSnext     = gampdf(DFcutS,chi0S,chi1S);
%             FSrho_plusSnext = mu_rhoS*(1-gamcdf(DFcutS,chi0S+1,chi1S));
%             FSrho_minusSnext= mu_rhoS*gamcdf(DFcutS,chi0S+1,chi1S);
            
            bCnext  = BCnext./KCnext;
            PiCnext = (1-eta)*Znext.*(nCnext).^eta  + pnext - ....
                deltaK + ((pnext-1).^2)./(2*phi) ;
            LCnext  = bCnext./PiCnext;
            DFcutC  = LCnext-deltaC -vCnext./PiCnext;
            FCnext  = gamcdf(DFcutC,chi0C,chi1C);
            FCrho_plusCnext  = mu_rhoC*(1-gamcdf(DFcutC,chi0C+1,chi1C));
            FCrho_minusCnext = mu_rhoC*gamcdf(DFcutC,chi0C+1,chi1C);
            
            phiICnext  = ( (phi/2)*((i_cnext - deltaK).^2) + phiK/2*(kCnext-1).^2  ).*KCnext ;
            
            DWL_Cnext = xiC*FCrho_minusCnext.* (PiCnext-(1-deltaK)*pnext).*KCnext ;
            Y_Cnext     = Znext.*(KCnext.^(1-eta)).*NCnext.^eta ;
            C_next      = Ynext + Y_Cnext ....
                - ICnext     - phiICnext ...
                - DWL_Cnext ;
            
            U1 = C.^(-gamma) ;
            U1next = C_next.^(-gamma) ;
            U2next = psi;
            
            SDF = beta * U1next./U1  ;
%             
%             if epsilon==0
%                 H_next=(BSnext.^alpha_eff) .* BCnext.^(1-alpha_eff);
%                 MRS_S_next = alpha_eff.*(U2next./U1next) .* H_next.^-gammaH .*(BCnext./BSnext).^(1-alpha_eff);
%                 MRS_C_next=  (1-alpha_eff).*(U2next./U1next) .* H_next.^-gammaH .*(BCnext./BSnext).^(-alpha_eff);
%             else
%                 H_next=( alpha_eff.*BSnext.^epsilon + (1-alpha_eff).*BCnext.^epsilon ).^(1/epsilon);
%                 MRS_S_next = alpha_eff.*(U2next./U1next) .* H_next.^-gammaH .*(H_next/BSnext).^(1-epsilon);
%                 MRS_C_next=  (1-alpha_eff).*(U2next./U1next) .* H_next.^-gammaH .*(H_next./BCnext).^(1-epsilon);                
%             end
            H_next = BCnext;
            MRS_C_next=  (U2next./U1next) .* H_next.^-gammaH;
            
            if params.use_kappa_fair
                FCrCnext=(1-xiC)*FCrho_minusCnext./LCnext;
                %                 FCrCnext=(1-xiC)*FCrho_minusCnext./LCnext;
                kappaC = prnext*(SDF.*(FCnext - FCrCnext));
            end
            
            % HH FOCs
            %FSrecov=(1-xiS)*FSrho_minusSnext.*(1-lSnext)./(LSnext.*(1-varrhonext));  % REV: new recovery for S-banks
            %fx(1) = qS - prnext*(SDF.*( (1-varrhonext).*( 1-FSnext+ FSnext*piB +(1-piB)*FSrecov) + varrhonext + MRS_S_next)); % REV: new HH FOC with runs
            fx(1) = qC - prnext*(SDF.*(1 + MRS_C_next));
            
            % shadow banks' FOCs
%             Lsptnext    = (1- lSnext./LSnext.*( deltaS*(1-xtnext) - xtnext.*vSnext./PiSnext) )./....
%                 ((1-lSnext.*(1-xtnext)).^2);
%             Lsptnext    = (1 - varrhonext)./(1-lSnext).^2;  % REV: new L script (app B.3.1)        
%             qSprimebS= -(1-piB)*prnext*(SDF.*( (1-xiS)*FSrho_minusSnext./(LSnext.*bSnext)....  % REV: new qSprime (app B.3)
%                 + Lsptnext.*(fSnext./bSnext).*(xiS*(1-varrhonext).*LSnext+(1-xiS)*(1-lSnext).*(deltaS + vSnext./PiSnext ))));
%             fx(3) = qS - kappaS + qSprimebS.*bSnext - lamSplus - prnext*(SDF.*( (1-FSnext).*( 1 - varrhonext + lSnext./LSnext.*vSnext./PiSnext)  ...
%                                                                                + FSrho_plusSnext.*lSnext./LSnext  ));  % REV: new FOC for qS
%             FStilde_next = FSrho_plusSnext.*(1-lSnext) - ...
%                 (1-FSnext).*((1-varrhonext).*LSnext - (1-lSnext).*vSnext./PiSnext) - FSnext.*(1-lSnext).*deltaS;  % REV: new expression for OmegaS
%             
%             fx(4) =  p + phiK*(kS-1) - (qS-kappaS)*bSnext - prnext*(SDF.*PiSnext.*FStilde_next)  ;
            
            % commercial banks' FOCs
            fx(2)        = qC - kappaC - lamCplus - prnext*(SDF.*(1-FCnext));
            FCtilde_next = FCrho_plusCnext - (1-FCnext).*(LCnext - vCnext./PiCnext) - FCnext*deltaC;
            fx(3) =  p + phiK*(kC-1) - (qC-kappaC)*bCnext - prnext*(SDF.*PiCnext.*FCtilde_next)  ;
            fx(4) = (1-theta)*p - bCnext - lamCminus;
            %fx(5) = thetaS*p - bSnext - lamSminus;
            
            % if mode==1, also compute marginal value functions
            V=[];
            if mode==1
                Vnext=zeros(obj.Vfct.Nof,1);
                VHnext= Pol_next(:,2);
                VH = (C^(1-gamma))./(1-gamma) + psi*(H^(1-gammaH))./(1-gammaH) + beta*prnext*VHnext;
                DC=envec.DC;
                DCnext=Pol_next(:,3);
                pCnext=Pol_next(:,4);
                pC = prnext*(SDF.*(DCnext+(1-FCnext).*pCnext));
                Vnext(1)=  p;
                Vnext(2)= VH;
                Vnext(3)= DC;
                Vnext(4)= pC;
                Vnext(5)= kC;
                V{1} = Vnext;
                V{2} = [];
            elseif mode==2
                % only during simulation
                
                % counterfactual risk free rate without liquidity services
                q = prnext*SDF;
                rf = 1/q-1;
                
                % conditional expected returns
                DCnext=Pol_next(:,3);
                pCnext=Pol_next(:,4);
                
                pC = prnext*(SDF.*(DCnext+(1-FCnext).*pCnext));
                exRC=prnext*(DCnext+(1-FCnext).*pCnext)/pC;  % not expected excess return
                exZ =prnext*PiCnext/p;
                
                
                if params.use_kappa_fair
                    kappa_fair=kappaC;
                else
                    FCrCnext=(1-xiC)*FCrho_minusCnext./LCnext;
                    kappa_fair=prnext*(SDF.*(FCnext - FCrCnext));
                end
                
                rets=struct('H',H,...
                    'C',C,...
                    'MRS_C',prnext*MRS_C_next,...
                    'rf',rf,...
                    'exRC',exRC,...
                    'exZ',exZ,...
                    'kappa_fair',kappa_fair,'theta',theta);
                
                
                V=rets;
            end
        end

    
        function [errmat,solmat,retmat]=calcEEError(obj,pointmat)
            % function to compute Euler equation error at points in state
            % space given by pointmat
            
            errmat=zeros(size(pointmat,1),obj.Pfct.Nof);
            solmat=zeros(size(errmat));
            retmat=zeros(size(pointmat,1),8);
            
            parfor i=1:size(errmat,1)
                point=pointmat(i,:);
                soltmp=obj.evaluatePol(point)';
                % transition
                [nextst,outstr]=obj.calcStateTransition(point,soltmp,0);
                % equations
                [fx,V]=obj.calcEquations(point(1),nextst,soltmp,outstr,2);                
                p = exp(soltmp(1));
                qC = exp(soltmp(2));
                normvec=[qC,p,qC,p];
                errmat(i,:)=fx'./normvec;
                solmat(i,:)=soltmp;
                retmat(i,:)=[V.H, V.C, V.MRS_C,...
                             V.rf, V.exRC, V.exZ, V.kappa_fair, V.theta];
            end
            
        end
        

        % simulate model
        function [simseries,varnames,errmat]=simulate(obj,NT,NTini,inistvec,simerror,shmat_in)
            if length(inistvec)~=obj.Vfct.SSGrid.Ndim
                error('inistvec must be vector of length SSGrid.Ndim');
            end
            
            NTtot=NT+NTini;
            simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD);
            
            if isempty(shmat_in)
                rng('default') % set seed
                rng(1); 
                shmat=rand(NTtot,1);
            else
                shmat=shmat_in;
            end
            point=inistvec;
            
            pointmat=zeros(NTtot,length(point)) ;
             for t=1:NTtot
                 
                try 
                    pointmat(t,:)=point;
                catch
                    disp(point);
                end
                exst=point(1);
                
                % next period's exog. state
                transprob=cumsum(obj.Exogenv.mtrans(exst,:));
                 exnext=find(transprob-shmat(t)>=0,1,'first');

                % transition to next period
                solvec=obj.evaluatePol(point)';
                valvec=obj.evaluateVal(point)';
                 
                [nextst,outstr]=obj.calcStateTransition(point,solvec,exnext);
                 
                addvec=model.DSGEModel.structToVec(outstr.addvars)';
                % write different categories of variables in one row
                simseries(t,:)=[point(1),outstr.exstvec',point(2:end),solvec,valvec,addvec];
                point=nextst;
            end
            
            simseries=simseries(NTini+1:end,:);
            varnames=[{'exst'}, obj.Ex_names, obj.En_names, ....
                  obj.Sol_names, obj.V_names, obj.Add_names];
             
            errmat=[];
            if simerror
                [errmat,~,retmat]=obj.calcEEError(pointmat);
                errmat=errmat(NTini+1:end,:);
                retmat=retmat(NTini+1:end,:);
                simseries=[simseries,retmat];
                varnames=[varnames,{'H','C','MRS_C','rf',....
                                  'exRC','exZ','kappa_fair','theta'}];       
                
               
            end
                        
        end
        
         function mobj=polIter(mobj,MAXIT,revisitFailed,printmode,avg_tol)
            gridSt=mobj.Vfct.SSGrid.Pointmat; % use points from BaseGrid here
            NPT=mobj.Vfct.SSGrid.Npt;
            exnpt=size(mobj.Exogenv.pts_perm,1);
            
            % initialize
            resmat=mobj.evaluatePol(gridSt)';
            resmat_prev=resmat;
            
            % split up matrix of points for better output
            gr_points = cell(exnpt,1);
            gr_index  = cell(exnpt,2);
            for i=1:exnpt
                grinlog=(gridSt(:,1)==i);
                grind=find(grinlog);
                gr_points{i}=gridSt(grinlog,:);
                gr_index{i,1}=grinlog;
                gr_index{i,2}=grind;
            end
            
            % value function
            VF=mobj.evaluateVal(gridSt)';
%            VF=mobj.evaluateVal(gridSt);
            VFnext=zeros(size(VF));
%             TF=zeros(size(VF,1),mobj.Tfct.Nof);
%             TFnext=TF;
                        
            % control flags
            iter=0;
                        
            disp(' ');
            disp('Starting main loop ...');
            disp(' ');
            while 1
                % counter
                iter=iter+1;
                
                % ===========================================
                % loop over state space
                % ===========================================
                
                % matrix for failed points
                failedPoints=[];
                % outer loop: all exogenous states
                for ei=1:exnpt
                    tmp_grid=gr_points{ei};
                    tmp_indlog=gr_index{ei,1};
                    tmp_index=gr_index{ei,2};
                    tmp_resmat=resmat(tmp_indlog,:);
                    tmp_resmat_prev=resmat_prev(tmp_indlog,:);
                    
                    disp(['State ',num2str(ei)]);

                    [tmp_resmat_new,tmp_VF,~,tmp_failed]=mobj.solvePointList(tmp_grid,tmp_resmat,tmp_resmat_prev,printmode,[]);
                    if revisitFailed
                        failedPoints=[failedPoints; tmp_index(tmp_failed)];
                    end
                                       
                    resmat_prev(tmp_indlog,:)=tmp_resmat;
                    resmat(tmp_indlog,:)=tmp_resmat_new;
                    VFnext(tmp_indlog,:)=tmp_VF;
                    %TFnext(tmp_indlog,:)=tmp_TF;
                end                          
                
                if ~isempty(failedPoints)
                    disp( '~~~~~~~~~~~~~~~~~~~');
                    disp(['Revisiting failed points: ',num2str(length(failedPoints)),' add. points ...']);
                    % try to solve at failed points
                    [new_resmat,new_VF,~,n_succ]=mobj.solvePointListFailed(gridSt,failedPoints,resmat,1,printmode,[]);
                    resmat(failedPoints,:)=new_resmat;
                    VFnext(failedPoints,:)=new_VF;
                    %TFnext(failedPoints,:)=new_TF;
                    disp(['Revisiting solved ',num2str(n_succ),' points.']);
                end                
                
                % approximate functions for next guess
                %VFnext=0.75*VFnext+0.25*VF;
                mobj=mobj.updateVfct(VFnext);
                %mobj=mobj.updateTfct(TFnext);
                
                  % convergence criterion (based on points in BaseGrid)
                val_range=2;
                VF_val = VF(:,val_range);
                VFnext_val=VFnext(:,val_range);
                [dist,wh]=max(abs(VF_val(:)-VFnext_val(:)));
                [mean_dist,col]=max(abs(mean(VF_val-VFnext_val)));
%                [distT,whT]=max(abs(TF(:)-TFnext(:)));                
                [wh_1,wh_2]=ind2sub(size(VFnext_val),wh);
 %               [whT_1,whT_2]=ind2sub(size(TFnext),whT);
                disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(mobj.V_names(val_range(1)-1+wh_2)), ...
                    ' at point ',num2str(wh_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:))]);
                disp(['-- Iteration: ',num2str(iter),', mean distance: ',num2str(mean_dist),' in ',char(mobj.V_names(val_range(1)-1+col))]);
%                 disp(['-- Iteration: ',num2str(iter),', max T distance: ',num2str(distT),' in col ',num2str(whT_2), ...
%                     ' at point ',num2str(whT_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(whT_1,:))]);
                disp(' ');
%                if dist<0.001 && mean_dist<0.0001 && distT<0.001
                if mean_dist<avg_tol
                    disp('Converged.');
                    break;
                elseif iter>=MAXIT
                    disp('Max.iter. exceeded.');
                    break;
                end
                
                % update guess
                VF=VFnext;
                %TF=TFnext;
            end
            
            % resulting policy functions
            mobj=mobj.updatePfct(resmat);       
         end
        
         
         function [simseries, varnames] = computeSimulationMoments(obj, simseries, varnames)

             % make HashMap with mapping of names to indices
             indexmap=java.util.HashMap;
             for i=1:length(varnames)
                 indexmap.put(varnames{i},i);
             end
             
             % list of indices
             loglist=model.HelperCollection.makeListFromNames(indexmap,{'qC'});
             multlist=model.HelperCollection.makeListFromNames(indexmap,{'lamC'});
             
             % conversion of log-values
             simseries(:,loglist)=exp(simseries(:,loglist));
             % conversion of multipliers
             simseries(:,multlist)=max(simseries(:,multlist),0).^(1/3);
             
             params=obj.Params;
             
             % ---------------------------------------------------------------------
             % state vars
             % ---------------------------------------------------------------------
             Y = simseries(:,indexmap.get('Y'));
             Z = simseries(:,indexmap.get('Z'));             
             Zbar = params.Zbar_disc*Z;
             varrho= simseries(:,indexmap.get('varrho'));
             
             % ---------------------------------------------------------------------
             % prices and interest rates
             % ---------------------------------------------------------------------
             p = simseries(:,indexmap.get('p'));
             qC = simseries(:,indexmap.get('qC'));
             rateC = 1./qC-1;
             EretZ=(Z(2:end)+p(2:end))./p(1:end-1);
             EEretZ_C=EretZ - rateC(1:end-1) -1;
             MRS_C = simseries(:,indexmap.get('MRS_C'));
             liqbenC =  MRS_C - params.kappaC;
             

             % ---------------------------------------------------------------------
             % C bank and S bank
             % ---------------------------------------------------------------------
             K = simseries(:,indexmap.get('K'));
             KC=K;
             BC = simseries(:,indexmap.get('bC')).*KC;
             deltaK=params.deltaK;
             deltaKbar=params.deltaKbar;
             Clev = BC./(KC.*p);
             Clevbk = BC./KC;
             nC = params.Nbar./KC;
             eta=params.eta;
             phi=params.phi;
             PiC = (1-eta)*Z.*nC.^eta + p - deltaK + ((p-1).^2)/(2*phi);
             LC = BC./(KC.*PiC);
             AC = BC;
             chi1C=params.sig_rhoC^2/params.mu_rhoC;
             chi0C=params.mu_rhoC/chi1C;
             
             % C bank
             kC = simseries(:,indexmap.get('kC'));
             vC = params.phiK/2*(kC.^2-1);
             FC = gamcdf(LC-params.deltaC-vC./PiC, chi0C, chi1C);
             rhoplus_C  = params.mu_rhoC*(1-gamcdf(LC-params.deltaC-vC./PiC, chi0C+1, chi1C))./(1-FC);
             FCrhominus_C = params.mu_rhoC*gamcdf(LC-params.deltaC-vC./PiC, chi0C+1, chi1C);
             DWL_C = params.xiC*FCrhominus_C.*PiC.*KC;
             lamC = simseries(:,indexmap.get('lamC'));
             rhominus_C=FCrhominus_C./FC;
             rhominus_C(FC==0)=0;
             recC = (1-params.xiC)*rhominus_C./LC;
             ERC = simseries(:,indexmap.get('exRC'));
             pC = simseries(:,indexmap.get('pC'));
             EERC = ERC - 1./qC;
             KCnext = KC .* kC;
             BCnext = simseries(:,indexmap.get('BCnext'));
             totvalC = pC + (qC-params.kappaC).*BCnext;
             rateC_eff = 1./(qC-params.kappaC) -1;
             waccC1 = (ERC-1).*pC./totvalC + rateC_eff.*(qC-params.kappaC).*BCnext./totvalC;
             waccC2 = (ERC-1).*(p.*KCnext-qC.*BCnext)./(p.*KCnext) + rateC_eff.*(qC.*BCnext)./(p.*KCnext);
             costperKC = p - (qC-params.kappaC).*BCnext./KC + params.phiK*(kC-1);

             
             assetsC= p.*KC;
            
              NC = nC.*KC ; 
  
              % Production
             YC =    Z.*(KC.^(1-eta)).*(NC.^eta) ;
             ishare = (p-1)/params.phi + deltaK;
             IC  = ishare.*KC ;
             I=  IC;
             irate = I./K;
            GDP = YC + Y ;
            FinShare = (YC)./GDP;
             krate = kC;
             pkrate = log(p(2:end).*K(2:end))-log(p(1:end-1).*K(1:end-1));
             iyrate = I./(YC);
             C=simseries(:,indexmap.get('C'));
             cyrate = C./GDP;
             logC = log(C);
             
             % ---------------------------------------------------------------------
             % HH and welfare
             % ---------------------------------------------------------------------
             Safe_sh = (BC)./(p.*K);
             DWL = DWL_C;
             
             % add to simseries
             simseries=[simseries(:,2:end),liqbenC,...
                                           rateC,EERC,waccC1,waccC2,costperKC,Clev,Clevbk,LC,FC,rhoplus_C,DWL_C,recC,assetsC,lamC,...
                                           I,irate,krate,iyrate,cyrate,logC,Safe_sh,YC, DWL, FinShare, GDP];
             
             simseries=[simseries(2:end,:),EretZ,EEretZ_C,  NC(2:end,:),....
                         AC(2:end,:), KC(2:end,:),pkrate ]; % dynamic variables  
                 
                 varnames_add={'liqbenC',...
                           'rateC','EERC','waccC1','waccC2','costperKC','Clev','Clevbk','LC','FC','rhoplus_C','DWL_C','recC','assetsC','lamC',...
                           'I','irate','krate','iyrate','cyrate','logC',...
                           'Safe_sh','YC','DWL','FinShare','GDP',...
                           'EretZ','EEretZ_C','NC','AC','KC','pkrate'};
                 
             
             varnames=[varnames(2:end), varnames_add];
             
         end
        
    end %of object methods
        
    
    %==============================================================================
    methods (Static)
        % static class-specific methods
        
      function [fx,stvals]=compStSt(x,params,printmode)
          
          % unpack parameters
          theta = params.theta;
          if isfield(params,'kappaC')
              kappaC     = params.kappaC ;
              kappaS     = params.kappaS ;
          else
              kappaC=params.kappa;
              kappaS=0;
          end
          nu    = params.nu;
          xiS   = params.xiS;
          xiC   = params.xiC;
          beta  = params.beta;
          psi   = params.psi;
          phi   = params.phi;
           
          gamma  = params.gamma;
          gammaH = params.gammaH;
          alpha = params.alpha;
          piB    = params.piB;
          if isfield(params,'epsilon')
              epsilon=params.epsilon;
          else
              epsilon=0;
          end
          
          deltaS = params.deltaS ;
          deltaC = params.deltaC ;
          deltaK = params.deltaK;
          deltaKbar = params.deltaKbar;
          
          eta    = params.eta;
          Nbar   = params.Nbar;
          
          Y       = params.muY;
          Z       = params.muZ;
          Zbar    = params.muZbar ; 
          varrhoS = params.varrhoS_ss ; 
          
          mu_rhoC  = params.mu_rhoC;
          mu_rhoS  = params.mu_rhoS;
          sig_rhoC = params.sig_rhoC;
          sig_rhoS = params.sig_rhoS;
          
          chi1C   = sig_rhoC^2/mu_rhoC;
          chi0C   = mu_rhoC/chi1C;
          chi1S   = sig_rhoS^2/mu_rhoS;
          chi0S   = mu_rhoS/chi1S;
                   
          % unpack variables
          if sum(~isreal(x))>0
                disp(x);
          end
          p  = exp(x(1));
          C  = exp(x(2));
          AC = exp(x(3));
%           K  = exp(x(4));

          
          phizero=false;
          if phi==0
                phi=1;
                phizero=true;
          end                    
          
          bC = (1-theta)*p;
          BC = AC;
          KC = BC/bC;      
          K = KC;
          
           NC = Nbar;
           nC = NC/KC ; 
         
          PiC = Z*(1-eta)*nC^eta + p-deltaK + (1/(2*phi))*(p-1)^2;
          LC  = bC/PiC;
           w  = eta*Z*nC^(eta-1);
          

          %  compute probs, rhos
          DFcutC = LC-deltaC ; 
          
          FC          = gamcdf(DFcutC, chi0C,chi1C);
          FCrho_plusC = mu_rhoC*(1-gamcdf(DFcutC, chi0C+1,chi1C));
          FCrho_minusC =mu_rhoC*gamcdf(DFcutC, chi0C+1,chi1C);
                   
          FC_spt = FCrho_plusC - (1-FC)*LC - FC*deltaC ;
                  
          
          % U1, U2, and H
          U1 = C.^(-gamma) ;
          U2 = psi;

          Hint=BC;
          H = 1/(1-gammaH)* (Hint^(1-gammaH));
          
          MRS_C=  (U2./U1) .* Hint.^-gammaH;
          
          % C-bank
           qC    = beta*(1+MRS_C);
           FCrC  = (1-xiC)*FCrho_minusC/LC;  
             rC  = FCrC/FC;
          
          
          % use fair kappa?
          if params.use_kappa_fair
                kappaC=beta*FC*(1-rC);
          end

          % solve for remaining quantities
          % also needs changin PiS and PiC respectively
          DC = FCrho_plusC*KC*PiC-(1-FC)*BC+(1-FC)*KC*((qC-kappaC)*bC-p);
          pC = (beta/(1-beta*(1-FC)))*DC;
          
          
          G = BC*(FC*(1-rC)-kappaC);
          M = beta ;
          
          lambC   = M*( MRS_C + FC ) - kappaC  ;
          W       = -Y - w + G + C + pC + qC*AC;

          VH      = (1/(1-M))*( (C^(1-gamma))./(1-gamma) + psi*H ) ;
          
            % consumption
            IC  = ((p-1)/phi  + deltaK)*KC ;
            i_c = IC/KC ;
            I   = IC; 
           phiI_C = (phi/2)*((i_c - deltaK)^2)*KC         ;
         
          DWL_C = xiC*FCrho_minusC* (PiC-(1-deltaK)*p)*KC ;
           Y_C = Z*(KC^(1-eta))*NC^eta ;
          DWL_C_total = xiC*FCrho_minusC* PiC *KC ;
           

          
          % equations
          if phizero
              fx(1)=p-1;
          else
              fx(1) = p - (qC-kappaC)*bC - beta*(PiC*FC_spt);
          end
          fx(2) = C - (Y  +   Y_C   ....
                      - IC     - phiI_C ...
                      - DWL_C);
          fx(3) = I - K + (1-deltaK)*(1-xiC*FCrho_minusC)*KC;
%          fx(4) = nC - (eta*Z/w)^(1/(1-eta));


          
          
          svals=[];
          
          % write output (and assign return structure)
          if printmode>=1
              rtC = rC*FC ;
             DepShare  = (BC)/(Y+Y_C );
             DepC      = (BC)/(C);
             % welfare consumption scale relative to benchVH
             scaleVH = ( C * H^(psi/(1-psi)) ) / ((1-M)*(1-gamma)*params.benchVH)^(1/((1-gamma)*(1-psi)));
             % not correct
             GDP=Y+Y_C;
              
              if printmode>1
                       

               disp('------- States -------');
                  disp(['K :',num2str(K)]);
                  disp(['bC :',num2str(bC)]);
                  
               disp('------- Capital -------');
                  disp(['GDP :',num2str(GDP)]);
                  disp(['FinDepGDP/GDP :',num2str((Y_C)/GDP)]);
                  disp(['p*K/GDP :',num2str(p*K/GDP)]);
                  disp(['I/K :',num2str(I/K)]);
                  disp('------- Labor -------');
                  disp(['Nall :',num2str(NC)]);
                  disp('------- Recovery -------');
                  disp(['rC :',num2str(rC)]);
                  disp('------- Leverage -------');
                  disp(['LC :',num2str(LC)]);
                  disp('------- Default Risk -------');
                  disp(['FC :',num2str(FC)]);
                  disp(['MRS_C :',num2str(MRS_C)]);
                  disp(['MRS_C - kappa :',num2str(MRS_C-params.kappaC)]);
                  disp(['SD(rhoC) :',num2str(PiC*params.sig_rhoC)]);
                  disp('------- Prices -------');
                  disp(['p :',num2str(p)]);
                  disp(['rateC :',num2str(100*(1/qC-1)) ' %']);
                  disp(['qC :',num2str(qC)]);
                      
                  disp('------- HH -------');
                  disp(['C / Y :',num2str(C/GDP)]);
                  disp(['C  :',num2str(C)]);
                  disp(['AC  :',num2str(BC)]);
                    disp(['Welfare :',num2str(VH)]);
%                   disp(['ScaleVH :',num2str(scaleVH)]);
                  disp('------- Gov -------');
                  disp(['G :',num2str(G)]);
                   
              
              end
 
              svals=struct;
              names={'K','KC','BC','bC','p','pC','qC','C','H','I',... 15
                  'G','Y','W','DC','M','rC','lambC','LC','DWL_total',...30
                  'U1','U2','MRS_C','FC','FCrho_plusC',...45
                  'FCrho_minusC','welfare','C_EV',...
                  'AC','w','PiC'}; %% 51
              
              
              vars=[K,KC,BC,bC,p,pC,qC,C,H,I,... 15
                  G,GDP,W,DC,M,rC,lambC,LC,DWL_C_total,...30
                  U1,U2,MRS_C,FC,FCrho_plusC,...45
                  FCrho_minusC,VH,scaleVH,...
                  AC,w,PiC];
              
              for i=1:length(names)
                  svals.(names{i})=vars(i);
              end
              
              
              Sol=struct( 'p',p,...
                  'qC',qC,...
                  'BCnext',BC,...
                  'lamC',lambC^(1/3));
              
              
              
              V=struct('p',p,...
                       'VH',VH,...
                       'DC',DC,...
                       'pC',pC,...
                       'kC',1);
                       
            
              Add =struct('Knext',K,...
                  'BCnext',BC,...
                  'kC',1,...
                  'H',H,...
                  'C',C,...
                  'DC',DC,...
                  'Z',Z,...
                  'Y',Y,....
                  'NC',NC);

              State=struct('K',K,...
                           'bC',bC);
              
              Sparam=struct('') ;
              %           'PhiB',PhiB,...
              %               'Mbar',MB,...
              %               'Zbar',ZA,...
              %               'dbar',dI,...
              %               'qbar',q);
              %
              statsout=struct([]) ; %'ndebt',nomdebt,...
              %               'norig',PhiB*bB,...
              %               'rB',rB,...
              %               'rD',rD,...
              %               'LTV',PhiB*aB/pkB,...
              %               'Lrate',lrate,...
              %               'PRrat',PRrat);
              
              stvals=struct('Sol',Sol,...
                  'V',V,...
                  'Add',Add,...
                  'State',State,...
                  'Sparam',Sparam,...
                  'statsout',statsout,...
                  'svals',svals);
          end
             
          
      end
                      
     
        function [solguessvec,Vguessvec,V_names]=assignGuess(stv)
            solguess=struct('p',log(stv.Sol.p),...
                'qC',log(stv.Sol.qC),...
                'BCnext',log(stv.Sol.BCnext),...
                'lamC',stv.Sol.lamC);
            
         
            solguessvec=model.DSGEModel.structToVec(solguess);
                        
            Vguess=struct('p',stv.Sol.p,...
                          'VH',stv.V.VH,...
                          'DC',stv.V.DC,...
                          'pC',stv.V.pC,...
                          'kC',stv.V.kC);
            
            Vguessvec=model.DSGEModel.structToVec(Vguess);
            V_names=fieldnames(Vguess);
            
        end

%      function [out1,out2]=compStSt_SP(x,params,printmode,compmode)
%           
%           % unpack parameters
%           thetaC = params.theta;
%           thetaS = 1-params.thetaS;
%           nu     = params.nu;
%           xiS    = params.xiS;
%           xiC    = params.xiC;
%           
%           beta   = params.beta;
%           eta    = params.eta ;
%           psi    = params.psi;
%           phi    = params.phi;
%           gamma  = params.gamma;
%           
%           alpha  = params.alpha;
%           deltaS = params.deltaS ;
%           deltaC = params.deltaC ;
%           deltaK = params.deltaK;
%           Nbar   = params.Nbar;
%           
%           Y      = params.muY;
%           Z      = params.muZ;
%           
%           mu_rhoC  = params.mu_rho(1);
%           mu_rhoS  = params.mu_rho(2);
%           sig_rhoC = params.sig_rho(1);
%           sig_rhoS = params.sig_rho(2);
%           
%           chi1C    = sig_rhoC^2/mu_rhoC;
%           chi0C    = mu_rhoC/chi1C;
%           chi1S    = sig_rhoS^2/mu_rhoS;
%           chi0S    = mu_rhoS/chi1S;
%                    
%           % unpack variables
%           LS = exp(x(1));
%           LC = exp(x(2)) ; 
%           KS = exp(x(3));
%           KC = exp(x(4));
%           NS = exp(x(5));
%           
%           
%           %  compute probs, rhos
%           DFcutC = LC-deltaC ; 
%           DFcutS = LS-deltaS ;
%            
%           FC           = gamcdf(DFcutC, chi0C,chi1C);
%           FCrho_minusC = mu_rhoC*gamcdf(DFcutC, chi0C+1,chi1C);
% %           fC           = gampdf(DFcutC, chi0C,chi1C);
%           fCtilde      = gampdf(DFcutC, chi0C+1,chi1C);
%           
%           FS           = gamcdf(DFcutS, chi0S,chi1S);
%           FSrho_minusS = mu_rhoS*gamcdf(DFcutS, chi0S+1,chi1S);
%           fS           = gampdf(DFcutS, chi0S,chi1S);
%           fStilde      = gampdf(DFcutS, chi0S+1,chi1S);
% 
%            M      = beta ;
%            IS     = KS*(deltaK + (1-deltaK)*xiS*FSrho_minusS) ;
%            IC     = KC*(deltaK + (1-deltaK)*xiC*FCrho_minusC) ;
%            I      = IS + IC;
%            
%            K      = KS+KC;         
%            if phi>0
%                pKS = 1+phi*(IS/KS-deltaK);
%                pKC = 1+phi*(IC/KC-deltaK);
%            else
%                pKS=1;
%                pKC=1;
%                phi=1;
%            end
%            NC     = Nbar - NS;
%            nC     = NC/KC;
%            nS     = NS/KS;
%            w      = eta*Z * nS^(eta-1) ;
%            
%            PiS    = Z*(1-eta)*nS^eta + pKS - deltaK + (1/(2*phi))*(pKS-1)^2;
%            PiC    = Z*(1-eta)*nC^eta + pKC - deltaK + (1/(2*phi))*(pKC-1)^2;
%            AS     = PiS*KS*LS;
%            AC     = PiC*KC*LC;
% 
%           Lambda_C = 1;
%           Lambda_S = (1-FS)^nu ;
%           alpha_eff = alpha * Lambda_S;
%            
%           DWL_S=xiS*FSrho_minusS*PiS*KS;
%           DWL_C=xiC*FCrho_minusC*PiC*KC;
%          
%               H = (1/(1-gammaH))*...
%                    (((BS^alpha_eff)* BC^(1-alpha_eff))^(1-gammaH)-1);
%            YC  =(1-xiC*FCrho_minusC)*Z*KC^(1-eta)*NC^eta;
%            YS  =(1-xiS*FSrho_minusS)*Z*KS^(1-eta)*NS^eta;
%            C   = Y + YC + YS - I  ....
%                    - phi/2*(IS/KS-deltaK)^2 *KS ....
%                    - phi/2*(IC/KC-deltaK)^2 *KC ;
%           
%            UC = C.^(-gamma) ;
%            UH = psi;
%            
%            VH      = (1/(1-M))*( (C^(1-gamma)-1)./(1-gamma) + psi*H )  ;
%          
%  
%           c1 = AC - (1-thetaC)*pKC*KC;
%           c2 = AS - (1-thetaS)*pKS*KS;
%           
%           
%           if strcmp(compmode,'fmin')
%                out1=-VH;
%                out2=[];
%           else
%               c(1)=c1;
%               c(2)=c2;
%               
%               out1=c;
%               out2=[];
%           end
%           
%           % write output (and assign return structure)
%           if printmode>=1
%                % output for graphs on determination of A^S and A^C
%                 GDP=Y+YS+YC;
%                 scaleVH = ( C * H^(psi/(1-psi)) ) / ((1-M)*(1-gamma)*params.benchVH)^(1/((1-gamma)*(1-psi)));
%                 
%               if printmode>1
%                        
%                   disp('------- Capital -------');
%                   disp(['K :',num2str(K)]);
%                   disp(['p*K/GDP :',num2str((pKS*KS+pKC*KC)/GDP)]);
%                   disp(['GDP :',num2str(GDP)]);
%                   disp(['I :',num2str(I)]);
%                   disp(['I/K :',num2str(I/K)]);
%                  % disp(['Z*K/GDP :',num2str(Z*K/(Z*K+Y))]);
%                   disp('------- Shadow Banks -------');
%                   disp(['KSsh :',num2str(KS/K)]);
%                   disp(['KS :',num2str(KS)]);
%                   disp(['NS :',num2str(NS)]);
%                   disp(['LS :',num2str(LS)]);
%                   disp(['FS :',num2str(FS)]);
%                     disp(['lambS :',num2str(c2)]);
% 
%                   disp('------- C Banks -------');
%                   disp(['KC :',num2str(KC)]);
%                   disp(['NC :',num2str(NC)]);
%                   disp(['LC :',num2str(LC)]);
%                   disp(['FC :',num2str(FC)]);
%                   disp(['lambC :',num2str(c1)]);
%                  disp('------- Liquidity -------');
%                   disp(['Lambda_S :',num2str(Lambda_S)]);
%                   disp('------- Prices -------');
%                   disp(['pKS :',num2str(pKS)]);
%                   disp(['pKC :',num2str(pKC)]);
%                   disp('------- HH -------');
%                   disp(['C :',num2str(C)]);
%                   disp(['C / Y :',num2str(C/GDP)]);
%                   disp(['H :',num2str(H)]);
%                   disp(['AC :',num2str(AC)]);
%                   disp(['AS :',num2str(AS)]);
%                   disp(['nS :',num2str(nS)]);
%                  disp(['Welfare :',num2str(VH)]);
%               end
%  
%               svals=struct;
%               names={'K','KS','KC','KSshare','NS','NC','AS','AC','Debt','pKS','pKC','C','H','I',... 15
%                    'Y','M','lambC','LS','LC','DWL_total',...30
%                     'UC','UH','LambdaS','FS','FC','fS',...45
%                   'FCrho_minusC','FSrho_minusS','welfare','C_EV'}; %% 51
%               
%               
%               vars=[K,KS,KC,KS/K,NS,NC,AS,AC,AS+AC,pKS,pKC,C,H,I,... 15
%                   GDP, M,c1 ,LS,LC,DWL_S+DWL_C,...30
%                   UC,UH,Lambda_S,FS,FC,fS,...45
%                   FCrho_minusC,FSrho_minusS,VH,scaleVH];
%               
%               for i=1:length(names)
%                   svals.(names{i})=vars(i);
%               end
%               
%               
%               
%               stvals=struct('svals',svals);
%               
%               out2=stvals;
%           end          
%       end
%          
             
        
        function [x,fx,exit]=tryOtherGuesses(fhand,gvec,options)
            % list of guesses
            
            %  Toggle Lagrange multipliers as follows:
            %
            %                        | only with muRG
            %   muRP   +   +   -   - |    -   -
            %   lamR   -   +   +   - |    +   -
            %   lamS   -   -   -   - |    -   -
            %   muSG   -   -   +   - |    +   -
            %   ---------------------+----------
            %   muRG   -   -   -   - |    +   + 
            
%            gindex={[9],[9],[5,9,10,11]};
%            gindex={[1,4],[1,4]};
            gindex={4};
            gvals={-.5};
            x=gvec;
            fx=ones(size(gvec));
            exit=0;
            
            for i=1:numel(gindex)
                newguess=gvec;
                newguess(gindex{i})=gvals{i};
                [x,fx,exit]=fsolve(fhand,newguess,options);
                if exit>0
                    break;
                end                                
            end
        end
        
    end % of static methods
    
    
end