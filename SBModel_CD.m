
classdef SBModel_CD < model.DSGEModel
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
        function obj=SBModel_CD(params,endogenv,exogenv,vfct,pfct,tfct)
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
            if isfield(params,'fullnu')
                fullnu=params.fullnu;
            else
                fullnu=0;
            end
            piB       = params.piB;
             % extract state variables
            exst=point(1);
            K=point(2);
            bS=point(3);
            bC=point(4);
            KSsh=point(5);
                        
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
            qS = exp(solvec(2));
            qC = exp(solvec(3));
            BSpol = exp(solvec(4));
            BCpol = exp(solvec(5));
            KSshpol = exp(solvec(6));
            nS     = exp(solvec(7)) ; 
            
            % small value functions for default cutoffs
            vS=phiK/2*(thisvals(8)^2-1);
            vC=phiK/2*(thisvals(9)^2-1);
             
            % this period
            KS = K*KSsh;
            KC = K-KS;
            BS = bS*KS;
            
            nH = nS*(Z/Zbar)^(1/(eta-1));
            
            PiS = (1-eta)*Z*nS^eta  + p - deltaK + ((p-1)^2)/(2*phi) ; 
             LS = bS/PiS ;
            PiH = (1-eta)*Zbar*nH^eta  + p *(1-deltaKbar) ; 
             xt = PiH./PiS;
             lS = varrho.*LS./xt ;
             NS = nS*KS*(1-lS) ;
             NH = nH*lS*KS ;
         DFcutS = ((1-varrho).*LS - (1-lS)*(vS/PiS + deltaS))./(1-lS) ;   % REV: new default cutoff
             FS = gamcdf(DFcutS,chi0S,chi1S);
           FSrho_minusS = mu_rhoS*gamcdf(DFcutS,chi0S+1,chi1S);
           FSrho_plusS  = mu_rhoS*(1-gamcdf(DFcutS,chi0S+1,chi1S));
           
             
             nC = nS ;
             NC = nC*KC ;
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
            IS  = ((p-1)/phi  + deltaK)*KS*(1-lS) ;  
            i_s = IS/((1-lS)*KS); 
            I   = (IC + IS) ; 
            % deadweight losses 
            DWL_S = xiS*FSrho_minusS* (1-lS)*(PiS-(1-deltaK)*p)*KS ; % REV: no run DWL on HH production
            DWL_C = xiC*FCrho_minusC* (PiC-(1-deltaK)*p)*KC ;
           % capital transition
            Kpol =   I  + (1-deltaK)*(1-xiC*FCrho_minusC)*KC ...
                        + (1-deltaK)*(1-xiS*FSrho_minusS)*KS*(1-lS)...
                        + (1-deltaKbar)*KS*lS ;
           KSpol = Kpol*KSshpol;
           KCpol = Kpol-KSpol;
           kS = KSpol/KS;
           kC = KCpol/KC;
            
           % consumption 
         phiI_S = (phi/2)*((i_s - deltaK)^2)*KS*(1-lS) + phiK/2*(kS-1)^2 *KS*(1-lS)  ; 
         phiI_C = (phi/2)*((i_c - deltaK)^2)*KC +        phiK/2*(kC-1)^2 *KC        ;
           Y_S = Z*((1-lS)*KS)^(1-eta)*NS^eta ;
           Y_C = Z*(KC^(1-eta))*NC^eta ;
           Y_H = Zbar*((lS*KS)^(1-eta))*NH^eta ;
           C    = Y  +  Y_S    +  Y_C    +  Y_H ....
                      - IC     - phiI_C ...
                      - IS     - phiI_S ....
                      - DWL_S  - DWL_C ;
                  % liquidity
                  %   Lambda_C = 1;
                  Lambda_S = ((1-varrho*fullnu)*(1-FS+piB*fullnu*FS))^nu ;
                  alpha_eff = alpha * Lambda_S;
                  
                  if epsilon==0
                      H = (BS.^alpha_eff).*(BC^(1-alpha_eff));
                  else
                      H = (alpha_eff*BS^epsilon + (1-alpha_eff)*BC^epsilon)^(1/epsilon) ;
                  end
                  
                  % dividends
                  DS = FSrho_plusS*KS*(1-lS)*PiS - ....  REV: only production from non-run capital
                      (1-FS)*(1-varrho)*BS + (1-FS)*((qS-kappaS)*BSpol -.... REV: only need to pay back deposits not already redeemed
                      p*KSpol - phiK/2*(kS-1)^2 *KS*(1-lS) );
                  DC = FCrho_plusC*KC*PiC - (1-FC)*BC + (1-FC)*((qC-kappaC)*BCpol - p*KCpol - phiK/2*(kC-1)^2 *KC );
                  
                  
                exnpt = obj.Exogenv.exnpt;

            if mode>0
    % simulation, mode contains number of next period's state
                cind=obj.Exogenv.pts_all(mode,end);
                Knext=Kpol;
				bSnext=BSpol/KSpol;
                bCnext=BCpol/KCpol;
                KSshnext=KSpol/Kpol;
            else
                % solution mode, compute next period's state variables for all
                % possible Markov states
                cind=obj.Exogenv.pts_all(:,end);
                Knext=Kpol*ones(exnpt,1);
				bSnext=BSpol/KSpol*ones(exnpt,1);
                bCnext=BCpol/KCpol*ones(exnpt,1);
                KSshnext=KSpol/Kpol*ones(exnpt,1);
            end
            
            nextst=[cind,Knext,bSnext,bCnext,KSshnext];            
            
            addvars=struct('Knext',Kpol,...
                'BSnext',BSpol,...
                'BCnext',BCpol,...
                'KSnext',KSpol,...
                'kS',kS,...
                'kC',kC,...
                'H',H,...
                'FS',FS,...
                'C',C,...
                'DS',DS,...
                'DC',DC,...
                'Z',Z,...
                'Y',Y,....
                'NH',NH,...
                'NC',NC,...
                'lS',lS,...
                'KS',KS);
            
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
            if isfield(params,'fullnu')
                fullnu=params.fullnu;
            else
                fullnu=0;
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
            qS   = exp(solvec(2));
            qC   = exp(solvec(3));
            nS   = exp(solvec(7)) ;
            lamC = solvec(8);
            lamS = solvec(9);
            
            lamCplus=max(0,lamC)^3;
            lamCminus=max(0,-lamC)^3;
            lamSplus=max(0,lamS)^3;
            lamSminus=max(0,-lamS)^3;
            
            
            % extract some state-dependent values
            envec = instr.addvars;
            Knext=envec.Knext;
            BSnext = envec.BSnext;
            BCnext = envec.BCnext;
            KSnext = envec.KSnext;
            C = envec.C;
            H = envec.H;
            NH = envec.NH;
            NC = envec.NC;
            lS = envec.lS;
            KS = envec.KS;
            KCnext = Knext-KSnext;
            kS = envec.kS;
            kC = envec.kC;
            if params.use_theta_var
                theta=instr.exstvec(4);
            end
            
            % compute expectation terms
            prnext = obj.Exogenv.mtrans(exst,:);
            Ynext  = obj.Exogenv.pts_perm(:,1) ;
            Znext  = obj.Exogenv.pts_perm(:,2) ;
            
            %Zbarnext=obj.Exogenv.pts_perm(:,3);
            Zbarnext   = Zbar_disc*Znext;
            varrhonext = obj.Exogenv.pts_perm(:,3);
            
            % projection evaluation
            Pol_next   = obj.evaluateVal(nextst)';
            pnext      = Pol_next(:,1);
            nSnext     = Pol_next(:,7);
            nHnext     = nSnext.*(Znext./Zbarnext).^(1/(eta-1));
            kSnext     = Pol_next(:,8);
            kCnext     = Pol_next(:,9);
            vSnext = phiK/2*(kSnext.^2 - 1);
            vCnext = phiK/2*(kCnext.^2 - 1);
            
            bSnext     = BSnext./KSnext;
            PiSnext    = (1-eta)*Znext.*(nSnext).^eta  + pnext - ....
                deltaK + ((pnext-1).^2)./(2*phi) ;
            LSnext     = bSnext./PiSnext ;
            
            PiHnext    = (1-eta)*Zbarnext.*(nHnext).^eta +  ....
                pnext *(1-deltaKbar) ;
            xtnext     = PiHnext./PiSnext;
            lSnext     = varrhonext.*LSnext./xtnext  ;
            
            
            nCnext = nSnext;
            NCnext = nCnext .* KCnext;
            NSnext     = nSnext .* (1-lSnext) .* KSnext;
            NHnext     = nHnext .* lSnext .* KSnext;
            
            % next period
            ICnext  = ((pnext-1)/phi  + deltaK).*KCnext ;
            i_cnext = ICnext/KCnext ;
            ISnext  = ((pnext-1)/phi  + deltaK).*KSnext.*(1-lSnext) ;
            i_snext = ISnext./((1-lSnext).*KSnext);
            
            DFcutS     = ((1-varrhonext).*LSnext-(1-lSnext).*(vSnext./PiSnext+deltaS))./(1-lSnext) ;  % REV: new default threshold
            FSnext     = gamcdf(DFcutS,chi0S,chi1S);
            fSnext     = gampdf(DFcutS,chi0S,chi1S);
            FSrho_plusSnext = mu_rhoS*(1-gamcdf(DFcutS,chi0S+1,chi1S));
            FSrho_minusSnext= mu_rhoS*gamcdf(DFcutS,chi0S+1,chi1S);
            
            bCnext  = BCnext./KCnext;
            PiCnext = (1-eta)*Znext.*(nCnext).^eta  + pnext - ....
                deltaK + ((pnext-1).^2)./(2*phi) ;
            LCnext  = bCnext./PiCnext;
            DFcutC  = LCnext-deltaC -vCnext./PiCnext;
            FCnext  = gamcdf(DFcutC,chi0C,chi1C);
            FCrho_plusCnext  = mu_rhoC*(1-gamcdf(DFcutC,chi0C+1,chi1C));
            FCrho_minusCnext = mu_rhoC*gamcdf(DFcutC,chi0C+1,chi1C);
            
            phiISnext  = ( (phi/2)*((i_snext - deltaK).^2) + phiK/2*(kSnext-1).^2 ).*KSnext.*(1-lSnext)  ;
            phiICnext  = ( (phi/2)*((i_cnext - deltaK).^2) + phiK/2*(kCnext-1).^2  ).*KCnext ;
            
            DWL_Snext = xiS*FSrho_minusSnext.* (1-lSnext).*(PiSnext-(1-deltaK)*pnext).*KSnext ;  % REV: no DWL on HH capital
            DWL_Cnext = xiC*FCrho_minusCnext.* (PiCnext-(1-deltaK)*pnext).*KCnext ;
            Y_Snext     = Znext.*(((1-lSnext).*KSnext).^(1-eta)).*NSnext.^eta ;
            Y_Cnext     = Znext.*(KCnext.^(1-eta)).*NCnext.^eta ;
            Y_Hnext     = Zbarnext.*((lSnext.*KSnext).^(1-eta)).*(NHnext.^eta);
            C_next      = Ynext + Y_Snext + Y_Cnext + Y_Hnext....
                - ICnext     - phiICnext ...
                - ISnext     - phiISnext ....
                - DWL_Snext  - DWL_Cnext ;
            
            % HH SDF
            Lambda_Snext = ((1-varrhonext*fullnu).*(1-FSnext+piB*fullnu*FSnext)).^nu ;
            alpha_eff = alpha * Lambda_Snext;
            
            U1 = C.^(-gamma) ;
            U1next = C_next.^(-gamma) ;
            U2next = psi;
            
            SDF = beta * U1next./U1  ;
            
            if epsilon==0
                H_next=(BSnext.^alpha_eff) .* BCnext.^(1-alpha_eff);
                MRS_S_next = alpha_eff.*(U2next./U1next) .* H_next.^-gammaH .*(BCnext./BSnext).^(1-alpha_eff);
                MRS_C_next=  (1-alpha_eff).*(U2next./U1next) .* H_next.^-gammaH .*(BCnext./BSnext).^(-alpha_eff);
            elseif epsilon==1
                H_next=( alpha_eff.*BSnext + (1-alpha_eff).*BCnext );
                MRS_S_next = alpha_eff.*(U2next./U1next) .* H_next.^-gammaH ;
                MRS_C_next=  (1-alpha_eff).*(U2next./U1next) .* H_next.^-gammaH ;              
            else
                H_next=( alpha_eff.*BSnext.^epsilon + (1-alpha_eff).*BCnext.^epsilon ).^(1/epsilon);
                MRS_S_next = alpha_eff.*(U2next./U1next) .* H_next.^-gammaH .*(H_next/BSnext).^(1-epsilon);
                MRS_C_next=  (1-alpha_eff).*(U2next./U1next) .* H_next.^-gammaH .*(H_next./BCnext).^(1-epsilon);                
            end
            
            if params.use_kappa_fair
                FCrCnext=(1-xiC)*FCrho_minusCnext./LCnext;
                %                 FCrCnext=(1-xiC)*FCrho_minusCnext./LCnext;
                kappaC = prnext*(SDF.*(FCnext - FCrCnext));
            end
            
            % HH FOCs
            FSrecov=(1-xiS)*FSrho_minusSnext.*(1-lSnext)./(LSnext.*(1-varrhonext));  % REV: new recovery for S-banks
            fx(1) = qS - prnext*(SDF.*( (1-varrhonext).*( 1-FSnext+ FSnext*piB +(1-piB)*FSrecov) + varrhonext + MRS_S_next)); % REV: new HH FOC with runs
            fx(2) = qC - prnext*(SDF.*(1 + MRS_C_next));
            
            % shadow banks' FOCs
            Lsptnext    = (1 - varrhonext)./(1-lSnext).^2;  % REV: new L script (app B.3.1)        
            qSprimebS= -(1-piB)*prnext*(SDF.*( (1-xiS)*FSrho_minusSnext./(LSnext.*bSnext)....  % REV: new qSprime (app B.3)
                + Lsptnext.*(fSnext./bSnext).*(xiS*(1-varrhonext).*LSnext+(1-xiS)*(1-lSnext).*(deltaS + vSnext./PiSnext ))));
            fx(3) = qS - kappaS + qSprimebS.*bSnext - lamSplus - prnext*(SDF.*( (1-FSnext).*( 1 - varrhonext + lSnext./LSnext.*vSnext./PiSnext)  ...
                                                                               + FSrho_plusSnext.*lSnext./LSnext  ));  % REV: new FOC for qS
            FStilde_next = FSrho_plusSnext.*(1-lSnext) - ...
                (1-FSnext).*((1-varrhonext).*LSnext - (1-lSnext).*vSnext./PiSnext) - FSnext.*(1-lSnext).*deltaS;  % REV: new expression for OmegaS
            
            fx(4) =  p + phiK*(kS-1) - (qS-kappaS)*bSnext - prnext*(SDF.*PiSnext.*FStilde_next)  ;
            
            % commercial banks' FOCs
            fx(5)        = qC - kappaC - lamCplus - prnext*(SDF.*(1-FCnext));
            FCtilde_next = FCrho_plusCnext - (1-FCnext).*(LCnext - vCnext./PiCnext) - FCnext*deltaC;
            fx(6) =  p + phiK*(kC-1) - (qC-kappaC)*bCnext - prnext*(SDF.*PiCnext.*FCtilde_next)  ;
            fx(7) = (1-theta)*p - bCnext - lamCminus;
            fx(8) = thetaS*p - bSnext - lamSminus;
            fx(9) = - Nbar + NC + nS*(1-lS)*KS + NH ;
            
            % if mode==1, also compute marginal value functions
            V=[];
            if mode==1
                Vnext=zeros(obj.Vfct.Nof,1);
                VHnext= Pol_next(:,2);
                VH = (C^(1-gamma))./(1-gamma) + psi*(H^(1-gammaH))./(1-gammaH) + beta*prnext*VHnext;
                DS=envec.DS;
                DC=envec.DC;
                DSnext=Pol_next(:,3);
                DCnext=Pol_next(:,4);
                pSnext=Pol_next(:,5);
                pCnext=Pol_next(:,6);
                pS = prnext*(SDF.*(DSnext+(1-FSnext).*pSnext));
                pC = prnext*(SDF.*(DCnext+(1-FCnext).*pCnext));
                Vnext(1)=  p;
                Vnext(2)= VH;
                Vnext(3)= DS;
                Vnext(4)= DC;
                Vnext(5)= pS;
                Vnext(6)= pC;
                Vnext(7)= nS;
                Vnext(8)= kS;
                Vnext(9)= kC;
                V{1} = Vnext;
                V{2} = [];
            elseif mode==2
                % only during simulation
                
                % counterfactual risk free rate without liquidity services
                q = prnext*SDF;
                rf = 1/q-1;
                
                % conditional expected returns
                DSnext=Pol_next(:,3);
                DCnext=Pol_next(:,4);
                pSnext=Pol_next(:,5);
                pCnext=Pol_next(:,6);
                
                pS = prnext*(SDF.*(DSnext+(1-FSnext).*pSnext));
                pC = prnext*(SDF.*(DCnext+(1-FCnext).*pCnext));
                exRS=prnext*(DSnext+(1-FSnext).*pSnext)/pS;  % note: only compute expected return here, 
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
                    'MRS_S',prnext*MRS_S_next,...
                    'rf',rf,...
                    'exRS',exRS,...
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
            retmat=zeros(size(pointmat,1),10);
            
            parfor i=1:size(errmat,1)
                point=pointmat(i,:);
                soltmp=obj.evaluatePol(point)';
                % transition
                [nextst,outstr]=obj.calcStateTransition(point,soltmp,0);
                % equations
                [fx,V]=obj.calcEquations(point(1),nextst,soltmp,outstr,2);                
                p = exp(soltmp(1));
                qS = exp(soltmp(2));
                qC = exp(soltmp(2));
                normvec=[qS,qC,qS,p,qC,p,p,p,p];
                errmat(i,:)=fx'./normvec;
                solmat(i,:)=soltmp;
                retmat(i,:)=[V.H, V.C, V.MRS_C, V.MRS_S,...
                             V.rf, V.exRS, V.exRC, V.exZ, V.kappa_fair, V.theta];
            end
            
        end
        

        % simulate model
        function [simseries,varnames,errmat,nextst]=simulate(obj,NT,NTini,inistvec,simerror,shmat_in,varargin)
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
            
            if nargin>6
                capdestshock=varargin{1};
            else
                capdestshock=[];
            end
            
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
                
                if ~isempty(capdestshock) && t==1
                    nextst(2)=nextst(2)*capdestshock;
                    nextst(3:4)=nextst(3:4)/capdestshock;
                end
                
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
                varnames=[varnames,{'H','C','MRS_C','MRS_S','rf','exRS',....
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
            VFnext=zeros(size(VF));
                        
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
                mobj=mobj.updateVfct(VFnext);
                
                  % convergence criterion (based on points in BaseGrid)
                val_range=2;
                VF_val = VF(:,val_range);
                VFnext_val=VFnext(:,val_range);
                [dist,wh]=max(abs(VF_val(:)-VFnext_val(:)));
                [mean_dist,col]=max(abs(mean(VF_val-VFnext_val)));
                [wh_1,wh_2]=ind2sub(size(VFnext_val),wh);
                disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(mobj.V_names(val_range(1)-1+wh_2)), ...
                    ' at point ',num2str(wh_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:))]);
                disp(['-- Iteration: ',num2str(iter),', mean distance: ',num2str(mean_dist),' in ',char(mobj.V_names(val_range(1)-1+col))]);
                disp(' ');
                if dist<avg_tol
                    disp('Converged.');
                    break;
                elseif iter>=MAXIT
                    disp('Max.iter. exceeded.');
                    break;
                end
                
                % update guess
                VF=VFnext;
            end
            
            % resulting policy functions
            mobj=mobj.updatePfct(resmat);       
         end
        
         
         function [simseries, varnames] = computeSimulationMoments(obj, simseries, varnames, varargin)

             % make HashMap with mapping of names to indices
             indexmap=java.util.HashMap;
             for i=1:length(varnames)
                 indexmap.put(varnames{i},i);
             end
             
             % list of indices
             loglist=model.HelperCollection.makeListFromNames(indexmap,{'qS','qC'});
             multlist=model.HelperCollection.makeListFromNames(indexmap,{'lamC','lamS'});
             
             % conversion of log-values
             simseries(:,loglist)=exp(simseries(:,loglist));
             % conversion of multipliers
             simseries(:,multlist)=max(simseries(:,multlist),0).^(1/3);
             
             if nargin>3
                firstrow=varargin{1};
                simseries=[firstrow; simseries];
             end
             
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
             qS = simseries(:,indexmap.get('qS'));
             qC = simseries(:,indexmap.get('qC'));
             rateS = 1./qS-1;
             rateC = 1./qC-1;
             EretZ=(Z(2:end)+p(2:end))./p(1:end-1);
             EEretZ_S=EretZ - rateS(1:end-1) -1;
             EEretZ_C=EretZ - rateC(1:end-1) -1;
             ratespr = 1./qS-1./qC;
             MRS_C = simseries(:,indexmap.get('MRS_C'));
             MRS_S = simseries(:,indexmap.get('MRS_S'));             
             MRSspr = MRS_C - MRS_S;
             liqbenC =  MRS_C - params.kappaC;
             

             % ---------------------------------------------------------------------
             % C bank and S bank
             % ---------------------------------------------------------------------
             K = simseries(:,indexmap.get('K'));
             KSsh = simseries(:,indexmap.get('KSsh'));
             lS = simseries(:,indexmap.get('lS'));
             KS=K.*KSsh ;
             KC=K-KS;
             BS = simseries(:,indexmap.get('bS')).*KS;
             BC = simseries(:,indexmap.get('bC')).*KC;
             BSsh = BS./(BS+BC);
             deltaK=params.deltaK;
             deltaKbar=params.deltaKbar;
             Slev = BS./(KS.*p);
             Slevbk = BS./KS;
             Clev = BC./(KC.*p);
             Clevbk = BC./KC;
             nS = simseries(:,indexmap.get('nS'));
             eta = params.eta;
             phi = params.phi;
             nC = nS;
             nH = nS .*(Zbar./Z).^(1/(1-eta));
             PiS = (1-eta)*Z.*nS.^eta + p - deltaK + ((p-1).^2)/(2*phi);
             PiC = (1-eta)*Z.*nC.^eta + p - deltaK + ((p-1).^2)/(2*phi);
             PiH = (1-eta)*Zbar.*nH.^eta + p*(1 - deltaKbar);
             LS = BS./(KS.*PiS);
             LC = BC./(KC.*PiC);
             AC = BC;
             AS = BS;
             chi1C=params.sig_rhoC^2/params.mu_rhoC;
             chi0C=params.mu_rhoC/chi1C;
             chi1S=params.sig_rhoS^2/params.mu_rhoS;
             chi0S=params.mu_rhoS/chi1S;
             
             % fire sales
             xt= PiH./PiS;
             ell= lS;
             
             % S bank
             kS = simseries(:,indexmap.get('kS'));
             vS = params.phiK/2*(kS.^2-1);
             DFcutS = ((1-varrho).*LS-(1-ell).*(vS./PiS+params.deltaS))./(1-ell);  % REV
             FS = gamcdf(DFcutS, chi0S, chi1S);
             rhoplus_S  = params.mu_rhoS*(1-gamcdf(DFcutS, chi0S+1, chi1S))./(1-FS);
             FSrhominus_S = params.mu_rhoS*gamcdf(DFcutS, chi0S+1, chi1S);
             DWL_S = params.xiS*FSrhominus_S.*PiS.*(1-ell).*KS;   % REV
             lamS_binds= (simseries(:,indexmap.get('lamS'))>0);  
             rhominus_S=FSrhominus_S./FS;
             rhominus_S(FS==0)=0;
             recS = (1-params.xiS)*rhominus_S./(LS.*(1-varrho));  % REV
             ERS = simseries(:,indexmap.get('exRS'));
             pS = simseries(:,indexmap.get('pS'));
             KSnext = KS .* kS;
             BSnext = simseries(:,indexmap.get('BSnext'));
             EERS = ERS - 1./qS;
             totvalS = pS + qS.*BSnext;
             waccS1 = (ERS-1).*pS./totvalS + rateS.*qS.*BSnext./totvalS;
             waccS2 = (ERS-1).*(p.*KSnext-qS.*BSnext)./(p.*KSnext) + rateS.*(qS.*BSnext)./(p.*KSnext);
             costperKS = p - qS.*BSnext./KS + params.phiK*(kS-1);
             
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
             capgrdiff = kC - kS;
             waccdiff1 = waccC1 - waccS1;
             waccdiff2 = waccC2 - waccS2;

             
             if ~isfield( params,'fullnu')
               params.fullnu= 1;  
             end
             
             LambdaS= ((1-varrho*params.fullnu).*(1-FS+FS*params.piB*params.fullnu)).^params.nu;           
             assetsC= p.*KC;
             assetsS= p.*KS;             
             
              NS = nS.*KS.*(1-lS) ;
              NC = nC.*KC ; 
              NH = nH.*lS.*KS ;
  
              % Production
             YH = Zbar.*((lS.*KS).^(1-eta)).*(NH.^eta) ;
             YS =    Z.*(((1-lS).*KS).^(1-eta)).*(NS.^eta) ;
             YC =    Z.*(KC.^(1-eta)).*(NC.^eta) ;
             YB = YS + YC;
            gYB = YB(2:end)./YB(1:end-1)-1;  
             ishare = (p-1)/params.phi + deltaK;
             IC  = ishare.*KC ;
             IS  = ishare.*(1-lS).*KS ;
             I=  IS + IC;
             irate = I./K;
            GDP = YH + YS + YC + Y ;
            FinShare = (YC + YS)./GDP;
             krate = (1-KSsh).*kC + KSsh.*kS;
             pkrate = log(p(2:end).*K(2:end))-log(p(1:end-1).*K(1:end-1));
             iyrate = I./(YS+YC);
             C=simseries(:,indexmap.get('C'));
             cyrate = C./GDP;
             logC = log(C);
             
             % Run-induced losses in ouput and excess depreciation
             DWL_run = (Z-Zbar).*((lS.*KS).^(1-eta)).*(NH.^eta) + (deltaKbar-deltaK)*lS.*KS;
             
             % ---------------------------------------------------------------------
             % HH and welfare
             % ---------------------------------------------------------------------
             Safe_sh = (BS+BC)./(p.*K);
             DWL = (DWL_S + DWL_C + DWL_run)./GDP;
             
             % add to simseries
             simseries=[simseries(:,2:end),rateS,EERS,waccS1,waccS2,costperKS,Slev,Slevbk,LS,FS,rhoplus_S,DWL_S,recS,assetsS,LambdaS,ratespr,MRSspr,liqbenC,...
                                           rateC,EERC,waccC1,waccC2,costperKC,Clev,Clevbk,LC,FC,rhoplus_C,DWL_C,recC,assetsC,lamC,lamS_binds,capgrdiff,...
                                           I,irate,krate,iyrate,cyrate,logC,Safe_sh,xt,lS,YH, YS, YC, YB, DWL, FinShare, BSsh, GDP, waccdiff1, waccdiff2];
             
             simseries=[simseries(2:end,:),EretZ,EEretZ_S,EEretZ_C, NS(2:end), NC(2:end,:),NH(2:end,:)....
                        AS(2:end), AC(2:end,:), KC(2:end,:),pkrate ,gYB]; % dynamic variables  
                 
                 varnames_add={'rateS','EERS','waccS1','waccS2','costperKS','Slev','Slevbk','LS','FS','rhoplus_S','DWL_S','recS','assetsS','LambdaS','ratespr','MRSspr','liqbenC',...
                           'rateC','EERC','waccC1','waccC2','costperKC','Clev','Clevbk','LC','FC','rhoplus_C','DWL_C','recC','assetsC','lamC','lamS_binds','capgrdiff',...
                           'I','irate','krate','iyrate','cyrate','logC',...
                           'Safe_sh','xt','lS','YH','YS','YC','YB','DWL','FinShare','BSsh','GDP','waccdiff1','waccdiff2',...
                           'EretZ','EEretZ_S','EEretZ_C','NS','NC','NH','AS','AC','KC','pkrate','gYB'};
                 
             
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
          bS = exp(x(2));
          C  = exp(x(3));
          AC = exp(x(4));
          AS = exp(x(5));
          nS = exp(x(6));
          K  = exp(x(7));
          muS = x(8);
          muC = x(9);
          if params.calib_eps
            epsilon=2/(1+exp(-x(10)))-1;
          end
          
          phizero=false;
          if phi==0
                phi=1;
                phizero=true;
          end                    
          
          bC = (1-theta)*p;
          BS = AS;
          BC = AC;
          KS = BS/bS;
          KC = BC/bC;        
          
          nH = nS*(Z/Zbar)^(1/(eta-1));
          PiS = (1-eta)*Z*nS^eta  + p - deltaK + ((p-1)^2)/(2*phi) ; 
           LS = bS/PiS ;
          PiH = (1-eta)*Zbar*nH^eta  + p *(1-deltaKbar) ; 
           xt = PiH./PiS;
           lS = varrhoS.*LS./xt ;
           NH = nH*lS*KS ;
           NS = nS*(1-lS)*KS;
           NC = Nbar - NS - NH ;
           nC = NC/KC ; 
         
          PiC = Z*(1-eta)*nC^eta + p-deltaK + (1/(2*phi))*(p-1)^2;
          LS  = bS/PiS;
          LC  = bC/PiC;
           w  = eta*Z*nS^(eta-1);
          
        
          muSplus=max(0,muS)^3;
          muSminus=max(0,-muS)^3;
          muCplus=max(0,muC)^3;
          muCminus=max(0,-muC)^3;
          

          %  compute probs, rhos
          DFcutC = real(LC-deltaC) ; 
          DFcutS = (1-varrhoS)*LS/(1-lS)-deltaS;  % REV: new default cutoff
          
          FC          = gamcdf(DFcutC, chi0C,chi1C);
          FCrho_plusC = mu_rhoC*(1-gamcdf(DFcutC, chi0C+1,chi1C));
          FCrho_minusC =mu_rhoC*gamcdf(DFcutC, chi0C+1,chi1C);
          
          FS          = gamcdf(DFcutS, chi0S,chi1S);
          FSrho_plusS = mu_rhoS*(1-gamcdf(DFcutS, chi0S+1,chi1S));
          FSrho_minusS= mu_rhoS*gamcdf(DFcutS, chi0S+1,chi1S);
          fS          = gampdf(DFcutS, chi0S,chi1S);
          
          FC_spt = FCrho_plusC - (1-FC)*LC - FC*deltaC ;
          FS_spt = FSrho_plusS*(1-lS) - (1-FS)*(1-varrhoS)*LS - FS*(1-lS)*deltaS ; % REV: new OmegaS
                  
          Lambda_S = (1-FS)^nu ;
          alpha_eff = alpha * Lambda_S;
          
          % U1, U2, and H
          U1 = C.^(-gamma) ;
          U2 = psi;
          if epsilon==0
              Hint=BS^alpha_eff * BC^(1-alpha_eff);
              H = 1/(1-gammaH)* (Hint^(1-gammaH));
              MRS_S = alpha_eff*(U2./U1) .* Hint.^-gammaH .*(BC./BS).^(1-alpha_eff);
              MRS_C=  (1-alpha_eff)*(U2./U1) .* Hint.^-gammaH .*(BC./BS).^(-alpha_eff);
          else
              Hint=( alpha_eff*BS^epsilon + (1-alpha_eff)*BC^epsilon )^(1/epsilon);
              H = 1/(1-gammaH)* (Hint^(1-gammaH));
              MRS_S = alpha_eff*(U2./U1) .* Hint.^-gammaH .*(Hint./BS).^(1-epsilon);
              MRS_C=  (1-alpha_eff)*(U2./U1) .* Hint.^-gammaH .*(Hint./BC).^(1-epsilon);              
          end                    
          
          % C-bank
           qC    = beta*(1+MRS_C);
           FCrC  = (1-xiC)*FCrho_minusC/LC;  
             rC  = FCrC/FC;
          
          % S-bank
            FSrS = (1-xiS)*(1-lS)*FSrho_minusS/(LS*(1-varrhoS));  % REV: new recovery
              rS = FSrS/FS;
              qS = beta*((1-varrhoS)*(1-FS + FS*piB+ (1-piB)*FSrS) + varrhoS + MRS_S);  % REV: new HH FOC for S-bank
          
          % use fair kappa?
          if params.use_kappa_fair
                kappaC=beta*FC*(1-rC);
          end

          % solve for remaining quantities
          % also needs changin PiS and PiC respectively
          DC = FCrho_plusC*KC*PiC-(1-FC)*BC+(1-FC)*KC*((qC-kappaC)*bC-p);
          pC = (beta/(1-beta*(1-FC)))*DC;
          
          DS=FSrho_plusS*KS*(1-lS)*PiS-(1-FS)*(1-varrhoS)*BS + (1-FS)*KS*((qS-kappaS)*bS-p); % REV: new S-bank dividend
          pS=(beta/(1-beta*(1-FS)))*DS;
          
          G = BC*(FC*(1-rC)-kappaC)+BS*piB*FS*(1-rS)-BS*kappaS;
          M = beta ;

          Lspt    = (1-varrhoS)/(1-lS)^2;  % REV: new L script
          qprimeS = - M*(1-piB)*( (1-xiS)*FSrho_minusS/(LS*bS) ....  % REV: new qSprime
                           + Lspt*(fS/bS)*((1-xiS)*(1-lS)*deltaS+xiS*(1-varrhoS)*LS)) ;
          
          lambC   = M*( MRS_C + FC ) - kappaC  ;
          W       = -Y - w + G + C + pC + pS + qS*AS + qC*AC;

          VH      = (1/(1-M))*( (C^(1-gamma))./(1-gamma) + psi*H ) ;
          
            % consumption
            IC  = ((p-1)/phi  + deltaK)*KC ;
            i_c = IC/KC ;
            IS  = ((p-1)/phi  + deltaK)*KS*(1-lS) ;  
            i_s = IS/((1-lS)*KS); 
            I   = (IC + IS) ; 
           phiI_S = (phi/2)*((i_s - deltaK)^2)*KS*(1-lS)  ; 
           phiI_C = (phi/2)*((i_c - deltaK)^2)*KC         ;
         
           DWL_S = xiS*FSrho_minusS* (1-lS)*(PiS-(1-deltaK)*p)*KS ;  % REV: no DWL on run capital
          DWL_C = xiC*FCrho_minusC* (PiC-(1-deltaK)*p)*KC ;
           Y_S = Z*(((1-lS)*KS)^(1-eta))*NS^eta ;
           Y_C = Z*(KC^(1-eta))*NC^eta ;
           Y_H = Zbar*((lS*KS)^(1-eta))*NH^eta ;
           DWL_S_total =  xiS*FSrho_minusS* (1-lS)*PiS*KS ;
          DWL_C_total = xiC*FCrho_minusC* PiC *KC ;
           

          
          % equations
          if phizero
              fx(1)=p-1;
          else
              fx(1) = p - (qC-kappaC)*bC - muCplus - beta*(PiC*FC_spt);
          end
          fx(2) = qS - kappaS + qprimeS*bS - beta*(1-FS)*(1-varrhoS) -  beta*lS*FSrho_plusS/LS;  % REV: new FOC for S banks
          fx(3) = C - (Y  +  Y_S    +  Y_C    +  Y_H ....
                      - IC     - phiI_C ...
                      - IS     - phiI_S ....
                      - DWL_S  - DWL_C);
          fx(4) = KC + KS - K;
          fx(5) = p - (qS-kappaS)*bS - muSplus - beta*(PiS*FS_spt);
          fx(6) = I - K + (1-deltaK)*(1-xiC*FCrho_minusC)*KC ...
                        + (1-deltaK)*(1-xiS*FSrho_minusS)*KS*(1-lS)...
                        + (1-deltaKbar)*KS*lS;
          fx(7) = nC - (eta*Z/w)^(1/(1-eta));
          % AS & AC <= 1
          fx(8)= KC  - muCminus;
          fx(9)= KS  - muSminus;
          % calibration mode
          if params.calib_eps
               fx(10)=params.betaScoef  - beta*((1-epsilon-gammaH)*(MRS_C/qC - MRS_S/qS )*alpha*(AS/Hint)^epsilon + ....
                        (1-epsilon).*(MRS_S/qS));
          end
          
          
          svals=[];
          
          % write output (and assign return structure)
          if printmode>=1
              rtS = rS*FS ;
              rtC = rC*FC ;
              % output for graphs on determination of A^S and A^C
              S_LHS = p-(qS-kappaS)*bS;
              S_RHS = beta*(PiS*FS_spt);
              C_LHS = p-(qC-kappaC)*bC;
              C_RHS = beta*(PiC*FC_spt) ;
              %      fS=fS*sig_rhoS; % actual density is std.normal density divided by sigma
               vS   = S_LHS - S_RHS ;
               vC   = C_LHS - C_RHS ;
             DepShare  = (BS + BC)/(Y+Y_S+Y_C + Y_H);
             DepC      = (BS + BC)/(C);
             % welfare consumption scale relative to benchVH
             scaleVH = ( C * H^(psi/(1-psi)) ) / ((1-M)*(1-gamma)*params.benchVH)^(1/((1-gamma)*(1-psi)));
             % not correct
             GDP=Y+Y_H+Y_C + Y_S;
              
              if printmode>1
                       

               disp('------- States -------');
                  disp(['K :',num2str(K)]);
                  disp(['bS :',num2str(bS)]);
                  disp(['bC :',num2str(bC)]);
                  disp(['KSsh :',num2str(KS/K)]);
                  
               disp('------- Capital -------');
                  disp(['GDP :',num2str(GDP)]);
                  disp(['FinDepGDP/GDP :',num2str((Y_S+Y_C)/GDP)]);
                  disp(['p*K/GDP :',num2str(p*K/GDP)]);
                  disp(['I/K :',num2str(I/K)]);
                   disp('------- Shadow Banks -------');
                  disp(['lS :',num2str(lS)]);
                  disp(['xt :',num2str(xt)]);
                  disp(['piB :',num2str(piB)]);
                  disp('------- Labor -------');
                  disp(['Nall :',num2str(NS + NC + NH)]);
                  disp('------- Recovery -------');
                  disp(['rS :',num2str(rS)]);
                  disp(['rC :',num2str(rC)]);
                  disp('------- Leverage -------');
                  disp(['LC :',num2str(LC)]);
                  disp(['Slev :',num2str(BS/(p*KS))]);
                  disp(['LS :',num2str(LS)]);
                  disp('------- Default Risk -------');
                  disp(['FS :',num2str(FS)]);
                  disp(['FC :',num2str(FC)]);
                  disp(['Lambda_S :',num2str(Lambda_S)]);
                  disp(['MRS_S :',num2str(MRS_S)]);
                  disp(['MRS_C :',num2str(MRS_C)]);
                  disp(['MRS_C - kappa :',num2str(MRS_C-params.kappaC)]);
                  disp(['SD(rhoC) :',num2str(PiC*params.sig_rhoC)]);
                  disp(['SD(rhoS) :',num2str(PiS*params.sig_rhoS)]);
                  disp('------- Prices -------');
                  disp(['p :',num2str(p)]);
                  disp(['qS :',num2str(qS)]);
                  disp(['rateC :',num2str(100*(1/qC-1)) ' %']);
                  disp(['qC :',num2str(qC)]);
                  disp(['spread % :',num2str( 100*(1./qS-1 -(1./qC-1)))]);
                  disp(['epsilon:',num2str(epsilon)]);
                      
                  disp('------- HH -------');
                  disp(['C / Y :',num2str(C/GDP)]);
                  disp(['C  :',num2str(C)]);
                  disp(['AC  :',num2str(BC)]);
                  disp(['AS  :',num2str(BS)]);
                  disp(['nS  :',num2str(nS)]);
                  disp(['Welfare :',num2str(VH)]);
                  disp('------- Gov -------');
                  disp(['G :',num2str(G)]);
                   
              
              end
 
              svals=struct;
              names={'K','KS','KC','KSshare','BS','BC','Debt','bS','bC','p','pS','pC','qS','qC','C','H','I',... 15
                  'vS','vC','G','Y','W','DS','DC','M','rC','lambC','LS','LC','DWL_total',...30
                  'rS','U1','U2','MRS_S','MRS_C','qprimeS','LambdaS','FS','FC','fS','FCrho_plusC',...45
                  'FCrho_minusC','FSrho_plusS','FSrho_minusS','welfare','C_EV',...
                  'S_LHS','S_RHS','C_LHS','C_RHS','Lspt','lS','xt','AS','AC','w','PiS','PiC','FS_spt'}; %% 51
              
              
              vars=[K,KS,KC,KS/K,BS,BC,BS+BC,bS,bC,p,pS,pC,qS,qC,C,H,I,... 15
                  vS,vC,G,GDP,W,DS,DC,M,rC,lambC ,LS,LC,DWL_S_total+DWL_C_total,...30
                  rS,U1,U2,MRS_S,MRS_C,qprimeS,Lambda_S,FS,FC,fS,FCrho_plusC,...45
                  FCrho_minusC,FSrho_plusS,FSrho_minusS,VH,scaleVH,...
                  S_LHS,S_RHS,C_LHS,C_RHS, Lspt,lS,xt, AS,AC,w,PiS,PiC,FS_spt];
              
              for i=1:length(names)
                  svals.(names{i})=vars(i);
              end
              
              
              Sol=struct( 'p',p,...
                  'qS',qS,...
                  'qC',qC,...
                  'BSnext',BS,...
                  'BCnext',BC,...
                  'KSshnext',KS/K,...
                  'nS',nS,....
                  'lamC',lambC^(1/3),...
                  'lamS',-0.3);
              
              
              
              V=struct('p',p,...
                       'VH',VH,...
                       'DS',DS,...
                       'DC',DC,...
                       'pS',pS,...
                       'pC',pC,...
                       'nS',nS,...
                       'kS',1,...
                       'kC',1);
              
            
              Add =struct('Knext',K,...
                  'BSnext',BS,...
                  'BCnext',BC,...
                  'KSnext',KS,...
                  'kS',1,...
                  'kC',1,...
                  'H',H,...
                  'FS',FS,...
                  'C',C,...
                  'DS',DS,...
                  'DC',DC,...
                  'Z',Z,...
                  'Y',Y,....
                  'NH',NH,...
                  'NC',NC,...
                  'lS',lS,...
                  'KS',KS);

              State=struct('K',K,...
                           'bS',bS,...
                           'bC',bC,...
                           'KSsh',KS/K);
              
              Sparam=struct('') ;

              statsout=struct([]); 
              
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
                'qS',log(stv.Sol.qS),...
                'qC',log(stv.Sol.qC),...
                'BSnext',log(stv.Sol.BSnext),...
                'BCnext',log(stv.Sol.BCnext),...
                'KSshnext',log(stv.Sol.KSshnext),...
                'nS',log(stv.Sol.nS),...
                'lamC',stv.Sol.lamC,...
                'lamS',stv.Sol.lamS);
            
         
            solguessvec=model.DSGEModel.structToVec(solguess);
                        
            Vguess=struct('p',stv.Sol.p,...
                          'VH',stv.V.VH,...
                          'DS',stv.V.DS,...
                          'DC',stv.V.DC,...
                          'pS',stv.V.pS,...
                          'pC',stv.V.pC,...
                          'nS',stv.V.nS,...
                          'kS',stv.V.kS,...
                          'kC',stv.V.kC);
            
            Vguessvec=model.DSGEModel.structToVec(Vguess);
            V_names=fieldnames(Vguess);
            
        end

        
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
            
            gindex={7,8};
            gvals={-.5,0.65};
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