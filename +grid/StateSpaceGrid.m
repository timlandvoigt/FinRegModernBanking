classdef StateSpaceGrid
   % abstract superclass
    
   properties (SetAccess=protected)
        StateBounds       
   end
   
   properties (Abstract, SetAccess=protected)
        Ndim
        Npt
        Pointmat 
        Type
   end
    
%--------------------------------------------------------------------------
% Regular methods (some checks for setter methods)
%--------------------------------------------------------------------------       
   methods
       function ssg=set.StateBounds(ssg,stBounds)
          nr=size(stBounds,1);
          if nr~=2 
              error('StateBounds must be a 2xNdim matrix');
          end
          ssg.StateBounds=stBounds;
      end
       
   end
%--------------------------------------------------------------------------
% End regular methods 
%--------------------------------------------------------------------------    


%--------------------------------------------------------------------------
% Static methods (miscellaenous functions for grids)
%--------------------------------------------------------------------------
   methods (Static)
       function nod=makeCombinations(x)
           
           dim=length(x);
           np=prod(x);
           nod=zeros(np,dim);
           
           % first dimension
           nod(:,1)=sort(repmat((1:x(1))',np/x(1),1));
           % remaining dims
           for d=2:dim
               % number of points divided by number of combinations from previous
               % dimensions
               rep_dim=np/prod(x(1:d));
               vec_elem=sort(repmat((1:x(d))',rep_dim,1));
               nod(:,d)=repmat(vec_elem,np/(rep_dim*x(d)),1);
           end          
       end
       
       
       function nod=makeCombinations_rev(x)
           
           dim=length(x);
           np=prod(x);
           nod=zeros(np,dim);
           
           % first dimension
           nod(:,dim)=sort(repmat((1:x(dim))',np/x(dim),1));
           % remaining dims
           for d=dim-1:-1:1
               % number of points divided by number of combinations from previous
               % dimensions
               rep_dim=np/prod(x(d:dim));
               vec_elem=sort(repmat((1:x(d))',rep_dim,1));
               nod(:,d)=repmat(vec_elem,np/(rep_dim*x(d)),1);
           end
       end

   end
%--------------------------------------------------------------------------
% End static methods 
%--------------------------------------------------------------------------
   
   
end