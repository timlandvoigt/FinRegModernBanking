classdef TensorGrid < grid.StateSpaceGrid
    
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        Ndim
        Npt
        Pointmat
        Type
        % tensor specific properties
        Unigrids    
        Dimvec
    end
    
    methods
        % constructor
        function ssg = TensorGrid(arg1,arg2)
            if nargin < 1
                error('Not enough input arguments.');
            end
            if nargin==1
                % if initialized by directly specifying univariate grids
                ssg=makeTensorGrid(ssg,arg1);
                ssg.StateBounds=zeros(2,ssg.Ndim);
                for d=1:ssg.Ndim
                    gridd=ssg.Unigrids{d};
                    ssg.StateBounds(:,d)=[gridd(1);gridd(end)];
                end
                ssg.Type='tensor (direct)';
            else
                % if initialized by specifying bounds and dimvec
                % construct equally space grid in that case
                ssg.StateBounds=arg1;
                dimvec=arg2;
                ndim=length(dimvec);
                unigrids=cell(ndim,1);
                for d=1:ndim
                    unigrids{d}=linspace(ssg.StateBounds(1,d),ssg.StateBounds(2,d),dimvec(d))';
                end
                ssg=makeTensorGrid(ssg,unigrids);
                ssg.Type='tensor (equ. sp.)';
            end
        end
        
        function obj = set.Unigrids(obj,gridarray)
            if (~iscell(gridarray) || length(size(gridarray))>2)
                error('grid array must be a (Ndim x 1) cell array');
            end
            if size(gridarray,2)>1
                gridarray=permute(gridarray,[2,1]);
            end
            if size(gridarray,2)>1
                error('grid array must be a (Ndim x 1) cell array');
            end
            obj.Unigrids=gridarray;
        end
               
        function grid=getUnigrid(ssg,dim)
           grid=ssg.Unigrids{dim}; 
        end        
        
        % make tensor grid from univariate grids
        function ssg = makeTensorGrid(ssg,unigrids)
            ssg.Unigrids = unigrids;
            ssg.Ndim=length(ssg.Unigrids);
            ssg.Dimvec=zeros(ssg.Ndim,1);
            for i=1:ssg.Ndim
                ssg.Dimvec(i)=length(ssg.Unigrids{i});
            end
            % generate point matrix
            ssg.Npt=prod(ssg.Dimvec);
            indmat=grid.StateSpaceGrid.makeCombinations_rev(ssg.Dimvec);
            ssg.Pointmat=zeros(ssg.Npt,ssg.Ndim);
            for d=1:ssg.Ndim
                gridd=ssg.Unigrids{d};
                ssg.Pointmat(:,d)=gridd(indmat(:,d));
            end
        end
        
        
    end
    
    
    
end