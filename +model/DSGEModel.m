classdef DSGEModel
% abstract superclass for model

    % common properties
    properties (Abstract, SetAccess=protected)
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
    
    methods (Abstract)
        calcStateTransition(obj,stvec,solvec,mode)
        calcEquations(obj,exst,nextst,solvec,envec,mode)
        calcEEError(obj,pointmat)
        polIter(obj,maxit)
    end
    
    methods
        % common constructor elements
        function obj=DSGEModel(params,endogenv,exogenv,vfct,pfct,tfct)
            if ~isa(pfct,'grid.ApproxFunction')
                error('pct must be a ApproxFunction.');
            end            
            obj.Pfct=pfct;
            if ~isa(vfct,'grid.ApproxFunction')
                error('vct must be a ApproxFunction.');
            end            
            obj.Vfct=vfct;      
            obj.NV=vfct.Nof;
            obj.V_names=reshape(endogenv.Vnames,1,obj.NV);
            obj.Tfct=tfct;
            if sum(isfield(endogenv,{'solnames','ennames','addnames'}))<3
                error('endogenv must contain cell arrays solnames, exnames, ennames and addnames as fields.');
            end            
            obj.NSOL=length(endogenv.solnames);
            obj.Sol_names=reshape(endogenv.solnames,1,obj.NSOL);
            if obj.NSOL~=pfct.Nof
                error('Inconsistent number of solution variables.');
            end
            obj.Sol_baseguess=endogenv.solbase;
            obj.NSTEN=length(endogenv.ennames);
            obj.En_names=reshape(endogenv.ennames,1,obj.NSTEN);
            if obj.NSTEN~=vfct.SSGrid.Ndim-1
                error('Inconsistent number of endogenous state variables.');
            end
            obj.NADD=length(endogenv.addnames);           
            obj.Add_names=reshape(endogenv.addnames,1,obj.NADD);            
            obj.Params=params;
            
            if sum(isfield(exogenv,{'exnames','mtrans','pts_perm','pts_all'}))<4
                error('exogenv must contain fields exnames, mtrans, pts_perm, and pts_all.');
            end
            obj.Exogenv=exogenv;
            obj.NSTEX=length(exogenv.exnames);
            obj.Ex_names=reshape(exogenv.exnames,1,obj.NSTEX);
            if obj.NSTEX~=size(exogenv.pts_perm,2)
                error('Inconsistent number of exogenous states');
            end
        end
        
        % evaluate policy functions
        function solvec=evaluatePol(obj,point)
            solvec=obj.Pfct.evaluateAt(point);
        end
        
        % evaluate 'value' functions
        function valvec=evaluateVal(obj,point)
            valvec=obj.Vfct.evaluateAt(point);
        end        
        
        % evaluate state transition function
        function stvec=evaluateTrans(obj,point)
            stvec=obj.Tfct.evaluateAt(point);
        end
        
        % solution at point based on abstract methods
        function [fx,V]=solveAtPoint(obj,solvec,point,mode,vals)
            if ~isreal(solvec)
               %disp(['Complex number at: ',num2str(point)]);
               solvec=real(solvec);
            end
            % transition
            [nextst,outstr]=obj.calcStateTransition(point,solvec,0,vals);
            % equations
            [fx,V]=obj.calcEquations(point(1),nextst,solvec,outstr,mode);
        end        
        
        function obj=updateVfct(obj,vals)
            obj.Vfct=obj.Vfct.fitTo(vals);
        end
        
        function obj=updatePfct(obj,vals)
            obj.Pfct=obj.Pfct.fitTo(vals);
        end
        
        function obj=updateTfct(obj,vals)
            obj.Tfct=obj.Tfct.fitTo(vals);
        end
        
        % simulate model 
        function [simseries,varnames,errmat]=simulate(obj,NT,NTini,inistvec,simerror)
            if length(inistvec)~=obj.Vfct.SSGrid.Ndim
                error('inistvec must be vector of length SSGrid.Ndim');
            end
            
            NTtot=NT+NTini;
            simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD);
        
            shmat=rand(NTtot,1);
            point=inistvec;
            
            pointmat=zeros(NTtot,length(point));
            
            for t=1:NTtot
               pointmat(t,:)=point; 
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
            
            errmat=[];
            if simerror
                errmat=obj.calcEEError(pointmat);
                errmat=errmat(NTini+1:end,:);
            end
            
            simseries=simseries(NTini+1:end,:);
            varnames=[{'exst'}, obj.Ex_names, obj.En_names, obj.Sol_names, obj.V_names, obj.Add_names];
            
        end
        
        % solve model at list of points
        function [resmat,VF,TF,failedPointsIndex]=solvePointList(mobj,tmp_grid,tmp_resmat,tmp_resmat_prev,print,logfname)
            
            gr_npt=size(tmp_grid,1);
            VF=zeros(gr_npt,mobj.Vfct.Nof);
            if ~isempty(mobj.Tfct)
                TF=zeros(gr_npt,mobj.Tfct.Nof);
            else
                TF=[];
            end
            resmat=tmp_resmat;
            failedPointsIndex=false(gr_npt,1);

            % solver options
            options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,'MaxIter',100,'MaxFunEvals',20000);            
          
            
            f=mobj.logopen(logfname);
            %fprintf(f,[repmat('.',1,gr_npt) '\n\n']);
            mobj.logclose(f);
            
            % Use "for" for debugging. Should be "parfor"
            parfor i=1:gr_npt
                ff=mobj.logopen(logfname);
                
                point=tmp_grid(i,:);
                
                succ='|';
                gvec=tmp_resmat(i,:)';
                % prepare function handle
                vals=mobj.evaluateVal(point);
                fh_solveAtPoint=@(sol)solveAtPoint(mobj,sol,point,0,vals);
                % call to fsolve
                [x,fx,exfl]=fsolve(fh_solveAtPoint,gvec,options);
                try
                    if exfl<1 && max(abs(gvec'-tmp_resmat_prev(i,:)))>0.01
                        [x,fx,exfl]=fsolve(fh_solveAtPoint,tmp_resmat_prev(i,:),options); 
                    end
                catch
                end
                if exfl<1
                    % retry with different guess
                    if print>1
                        fprintf(ff,['Try other guesses at point ',num2str(i),':',num2str(point),'\n']);
                    end
                    [xret,fx,exflret]=mobj.tryOtherGuesses(fh_solveAtPoint,gvec,options);
                    if exflret>0
                        xr=xret;
                        succ='t';
                    else
                        if print>0
                            fprintf(ff,['Return code: ',num2str(exfl),' at point ' ,num2str(i),':',num2str(point),'\n']);
                        end
                        xr=gvec;
                        succ='x';
                        failedPointsIndex(i)=true;
                    end
                else
                    xr=x;
                end
                
                % record result
                resmat(i,:)=xr;
                
                % also compute VF at solution
                [~,V]=solveAtPoint(mobj,xr,point,1,vals);
                
                VF(i,:)=V{1};
                if ~isempty(mobj.Tfct)
                    TF(i,:)=V{2};  
                end
                %fprintf(ff,'\b%s\n',succ);
                
                mobj.logclose(ff);
            end  
            
        end

        % solve model at list of points
        function [resmat_failed,VF_failed,TF_failed,n_succ]=solvePointListFailed(mobj,grid,index_failed,resmat,npasses,print,logfname)
            
            gr_npt=size(grid,1);
            VF=zeros(gr_npt,mobj.Vfct.Nof);
            if ~isempty(mobj.Tfct)
                TF=zeros(gr_npt,mobj.Tfct.Nof);
                TF_failed=zeros(length(index_failed),size(TF,2));
                updTF=true;
            else
                TF_failed=[];
                updTF=false;
            end
            resmat_failed=zeros(length(index_failed),size(resmat,2));
            VF_failed=zeros(length(index_failed),size(VF,2));

            idx_succ = zeros(length(index_failed));
            
            % solver options
            options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,'MaxIter',100,'MaxFunEvals',20000);            
          
            % keep track of failed indices
            failed_idx_list=1:length(index_failed);
                        
            for p=1:npasses
                
                new_index_failed=false(size(index_failed));              
                tmp_index_worked = setdiff(1:gr_npt,index_failed);
                tmp_grid_worked=grid(tmp_index_worked,:);
                
                f=mobj.logopen(logfname);
                fprintf(f,[repmat('.',1,length(index_failed)) '\n\n']);
                mobj.logclose(f);
                
                this_rmfailed=zeros(length(index_failed),size(resmat,2));
                % first loop over all failed points and find closest successful
                % point
                parfor i=1:length(index_failed)
                    ff=mobj.logopen(logfname);                  
                    
                    this_index=index_failed(i);
                    point=grid(this_index,:);
                    % Find closest state
                    dist=sum((tmp_grid_worked-repmat(point,length(tmp_grid_worked),1)).^2,2);
                    [~,distidx]=sort(dist);
					ntr=1+2*p;
                    gveclist = resmat(tmp_index_worked(distidx(1:ntr)),:);
                    
                    % prepare function handle
                    vals=mobj.evaluateVal(point);
                    fh_solveAtPoint=@(sol)solveAtPoint(mobj,sol,point,0,vals);
                    for g=1:ntr
                        gvec=gveclist(g,:);
                        % call to fsolve
                        try
                            [x,fx,exfl]=fsolve(fh_solveAtPoint,gvec,options);
                            if exfl>0
                                break;
                            end
                        catch
                            exfl=0;
                        end
                    end
                    if exfl>0
                        xr=x;
                        succ='|';
                        idx_succ(i)=1;
                    else
                        succ='x';
                        new_index_failed(i)=true;
                        xr=gveclist(1,:);
                    end
                    
                    % record result
                    this_rmfailed(i,:)=xr;
                    
                    % also compute VF at solution
                    [~,V]=solveAtPoint(mobj,xr,point,1,vals);
                    
                    VF_failed(i,:)=V{1};
                    if updTF
                        TF_failed(i,:)=V{2};
                    end
                    %fprintf(ff,'\b%s\n',succ);
                    
                    mobj.logclose(ff);
                end
                % write into return object
                resmat_failed(failed_idx_list,:)=this_rmfailed;
                failed_idx_list=failed_idx_list(new_index_failed);
                
                % update total resmat for gueses
                if sum(new_index_failed)==0 || sum(new_index_failed)==length(index_failed)
                    break;
                end
                resmat(index_failed(~new_index_failed),:)=resmat_failed(~new_index_failed,:);
                index_failed=index_failed(new_index_failed);
            
            end
            
            n_succ = sum(idx_succ);
        end
        
                  
    end

%--------------------------------------------------------------------------
% Static methods (miscellaenous functions for grids)
%--------------------------------------------------------------------------
   methods (Static)
       
       function f=logopen(logfname,perm)
          if nargin<2
              perm='a';
          end
          if isempty(logfname)
                f=1;
            else
                f=fopen(logfname,perm);
          end 
       end
       
       function logclose(f)
          if f>1
              fclose(f);
          end
       end
       
       function st=vecToStruct(vec,names)
           st=num2cell(vec);
           st=cell2struct(st,names,1);
       end
       
       function vec=structToVec(st)
            vec=struct2cell(st);
            vec=cell2mat(vec);
       end
       
       function [P_Rouw, z_Rouw] = rouwen(rho_Rouw, mu_uncond, sig_uncond, S_Rouw, n_R)
           %ROUWEN   Rouwenhorst's method (1995) to approximate an AR(1) process using
           %   a  finite state Markov process.
           %
           %   For details, see Rouwenhorst, G., 1995: Asset pricing  implications of
           %   equilibrium business cycle models, in Thomas Cooley (ed.), Frontiers of
           %   Business Cycle Research, Princeton University Press, Princeton, NJ.
           %
           %   Suppose we need to approximate the following AR(1) process:
           %
           %                   y'=rho_Rouw*y+e
           %
           %   where abs(rho_Rouw)<1, sig_uncond=std(e)/sqrt(1-rho_Rouw^2) and
           %   mu_uncond denotes E(y), the unconditional mean of y. Let n_R be the
           %   number of grid points. n_R must be a positive integer greater than one.
           %
           %   [P_Rouw, z_Rouw] = rouwen(rho_Rouw, mu_uncond, sig_uncond, n_R) returns
           %   the discrete state space of n_R grid points for y, z_Rouw, and
           %   the centrosymmetric transition matrix P_Rouw. Note that
           %
           %       1. z_Rouw is a column vector of n_R real numbers.
           %       2. The (i,j)-th element of P_Rouw is the conditional probability
           %          Prob(y'=z_Rouw(i)|y=z_Rouw(j)), i.e.
           %
           %                 P_Rouw(i,j)=Prob(y'=z_Rouw(i)|y=z_Rouw(j))
           %
           %           where z_i is the i-th element of vector z_Rouw. Therefore
           %
           %           P_Rouw(1,j)+P_Rouw(2,j)+ ... +P_Rouw(n,j)=1 for all j.
           %
           %   See also HITM_Z and HITM_S on how to simulate a Markov processes using
           %   a transition matrix and the grids.
           %
           %   Damba Lkhagvasuren, June 2005
           
           % CHECK IF abs(rho)<=1
           if abs(rho_Rouw)>1
               error('The persistence parameter, rho, must be less than one in absolute value.');
           end
           
           % CHECK IF n_R IS AN INTEGER GREATER THAN ONE.
           if n_R <1.50001 %| mod(n_R,1)~=0
               error('For the method to work, the number of grid points (n_R) must be an integer greater than one.');
           end
           
           % CHECK IF n_R IS AN INTEGER.
           if mod(n_R,1)~=0
               warning('the number of the grid points passed to ROUWEN is not an integer. The method rounded n_R to its nearest integer.')
               n_R=round(n_R);
               disp('n_R=');
               disp(n_R);
           end
           
           skewfun=@(p)(S_Rouw - (2*p-(1+rho_Rouw))/sqrt((n_R-1)*(1-p)*(p-rho_Rouw)));
           options=optimset('Display','none');
           [p,~,exit]=fsolve(skewfun,(rho_Rouw+1)/2,options);
           if exit<=0
              error('Could not solve for p given skewness'); 
           else
              q=rho_Rouw + 1 - p;
              if q<0 || q>1
                  error('Level of skewness not attainable given number of nodes'); 
              end
           end

           % GRIDS
           step_R = sig_uncond*sqrt((n_R-1)*(2-p-q)^2/(4*(1-p)*(1-q)));
           M=mu_uncond-step_R*(q-p)/(2-p-q);
           z_Rouw=[-1:2/(n_R-1):1]';
           z_Rouw=M+step_R*z_Rouw;
           
           % CONSTRUCTION OF THE TRANSITION PROBABILITY MATRIX           
           P_Rouw=[ p  (1-p);
               (1-q) q];
           
           for i_R=2:n_R-1
               a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
               a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
               a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
               a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
               P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
               P_Rouw(2:i_R, :) = P_Rouw(2:i_R, :)/2;
           end
           
           P_Rouw=P_Rouw';
           
           for i_R = 1:n_R
               P_Rouw(:,i_R) = P_Rouw(:,i_R)/sum(P_Rouw(:,i_R));
           end
       end
           
       
       function des = hitm_s(PPP,a1)
           %HITM_S simulates a time series for states (not actual values) of a finite
           %   state Markov chain using its transition matrix, PPP, and a column
           %   vector of N random numbers drawn from the standard uniform distribution
           %   on the interval(0,1), a1.
           %
           %   The length of the simulated time series is given by the size of a1,
           %   i.e. if a1 is an Nx1 column vector of the random numbers, the size of
           %   des will be Nx1.
           %
           %                   des = hitm_s(PPP,a1)
           %
           %   where
           %
           %       PPP - the transition matrix,
           %       a1 - a column vector of random numbers drawn from the standard
           %            uniform distribution on the interval(0,1), and
           %       des - the simulated time series for the states. If PPP is an MxM
           %             matrix, each period des will take one of the following M
           %             integers, {1,2,...,M}.
           %
           %   Note that the method assumes that (i,j)-th element of PPP is the
           %   conditional probability Prob(state(t+1)=i|state(t)=j), i.e.
           %
           %                 PPP(i,j)=Prob(state(t+1)=i|state(t)=j)
           %
           %   where i and j denote the numbers of the states. Therefore
           %
           %           PPP(1,j)+PPP(2,j)+ ... +PPP(n,j)=1, for all j.
           %
           %   See also HITM_Z for simulating actual values.
           %   Damba Lkhagvasuren, June 2005
           
           N=size(a1,1);
           znum=size(PPP,1);
           
           A=tril(ones(znum))*PPP;
           A(znum,:)=2;
           
           des=zeros(N,1);
           
           ainit=randperm(znum);
           des(1,1)=ainit(1,1);
           destemp=des(1,1);
           
           for c_ount=2:N;
               
               if a1(c_ount,1)<=A(1,destemp);
                   des(c_ount,1)=1;
               end ;
               
               for i=1:znum-1;
                   if A(i,destemp)<a1(c_ount,1)
                       if A(i+1,destemp)>=a1(c_ount,1);
                           des(c_ount,1)=i+1;
                       end
                   end
               end
               destemp=des(c_ount,1);
               
           end;
       end
       
       function [Z,Zprob] = tauchenhussey(N,mu,rho,sigma,baseSigma)
           % Function tauchenhussey
           %
           % Purpose:    Finds a Markov chain whose sample paths
           %             approximate those of the AR(1) process
           %                 z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
           %             where eps are normal with stddev sigma
           %
           % Format:     {Z, Zprob} = TauchenHussey(N,mu,rho,sigma,m)
           %
           % Input:      N         scalar, number of nodes for Z
           %             mu        scalar, unconditional mean of process
           %             rho       scalar
           %             sigma     scalar, std. dev. of epsilons
           %             baseSigma scalar, std. dev. used to calculate Gaussian
           %                       quadrature weights and nodes, i.e. to build the
           %                       grid. I recommend that you use baseSigma = w*sigma +
           %                       (1-w)*sigmaZ where sigmaZ = sigma/sqrt(1-rho^2),
           %                       and w = 0.5 + rho/4. Tauchen & Hussey recommend
           %                       baseSigma = sigma, and also mention baseSigma = sigmaZ.
           %
           % Output:     Z       N*1 vector, nodes for Z
           %             Zprob   N*N matrix, transition probabilities
           %
           %     Martin Floden, Stockholm School of Economics
           %     January 2007 (updated August 2007)
           %
           %     This procedure is an implementation of Tauchen and Hussey's
           %     algorithm, Econometrica (1991, Vol. 59(2), pp. 371-396)
           
           Z     = zeros(N,1);
           Zprob = zeros(N,N);
           
           [Z,w] = gaussnorm(N,mu,baseSigma^2);   % See note 1 below
           
           
           for i = 1:N
               for j = 1:N
                   EZprime    = (1-rho)*mu + rho*Z(i);
                   Zprob(i,j) = w(j) * norm_pdf(Z(j),EZprime,sigma^2) / norm_pdf(Z(j),mu,baseSigma^2);
               end
           end
           
           for i = 1:N
               Zprob(i,:) = Zprob(i,:) / sum(Zprob(i,:),2);
           end
           
           
           
           function c = norm_pdf(x,mu,s2)
               c = 1/sqrt(2*pi*s2) * exp(-(x-mu)^2/2/s2);
           end
           
           function [x,w] = gaussnorm(n,mu,s2)
               % Find Gaussian nodes and weights for the normal distribution
               % n  = # nodes
               % mu = mean
               % s2 = variance
               
               [x0,w0] = gausshermite(n);
               x = x0*sqrt(2*s2) + mu;
               w = w0 / sqrt(pi);
           end
           
           function [x,w] = gausshermite(n)
               % Gauss Hermite nodes and weights following "Numerical Recipes for C"
               
               MAXIT = 10;
               EPS   = 3e-14;
               PIM4  = 0.7511255444649425;
               
               x = zeros(n,1);
               w = zeros(n,1);
               
               m = floor(n+1)/2;
               for i=1:m
                   if i == 1
                       z = sqrt((2*n+1)-1.85575*(2*n+1)^(-0.16667));
                   elseif i == 2
                       z = z - 1.14*(n^0.426)/z;
                   elseif i == 3
                       z = 1.86*z - 0.86*x(1);
                   elseif i == 4
                       z = 1.91*z - 0.91*x(2);
                   else
                       z = 2*z - x(i-2);
                   end
                   
                   for iter = 1:MAXIT
                       p1 = PIM4;
                       p2 = 0;
                       for j=1:n
                           p3 = p2;
                           p2 = p1;
                           p1 = z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3;
                       end
                       pp = sqrt(2*n)*p2;
                       z1 = z;
                       z = z1 - p1/pp;
                       if abs(z-z1) <= EPS, break, end
                   end
                   if iter>MAXIT, error('too many iterations'), end
                   x(i)     = z;
                   x(n+1-i) = -z;
                   w(i)     = 2/pp/pp;
                   w(n+1-i) = w(i);
               end
               x(:) = x(end:-1:1);
           end
           
           % Note 1: If you have Miranda and Fackler's CompEcon toolbox you can use
           % their qnwnorm function to obtain quadrature nodes and weights for the
           % normal function: [Z,w] = qnwnorm(N,mu,baseSigma^2);
           % Compecon is available at http://www4.ncsu.edu/~pfackler/compecon/
           % Otherwise, use gaussnorm as here.
           
       end
   

       
       function [stlist,transmat,indlist]=makeTotalTransition(pointList,probList)
           
           N=1;
           NI=[];
           ndim=length(pointList);
           
           % total size of state space
           for i=1:ndim
               ni=length(pointList{i});
               N=N*ni;
               NI=[NI,ni];
           end
           
           indlist=grid.StateSpaceGrid.makeCombinations(NI);
           
           stlist=zeros(N,ndim);
           for d=1:ndim
               stlist(:,d)=pointList{d}(indlist(:,d));
           end
           
           transmat=zeros(N);
           for i=1:N
               dind=indlist(i,:);
               for j=1:N
                   tind=indlist(j,:);
                   tp=1;
                   for d=1:ndim
                       tp=tp*probList{d}(dind(d),tind(d));
                   end
                   transmat(i,j)=tp;
               end
           end                      
       end
       
       
   end
%--------------------------------------------------------------------------
% End static methods 
%--------------------------------------------------------------------------
    
    
    
end