classdef ApproxFunction
% abstract superclass for approximating function   
    
    % common properties
    properties (Abstract, SetAccess=protected)
        SSGrid % state space grid
        Nof % # functions
        Vals % matrix with interpolated/approximated data points
    end
    
    % common methods to be implemented
    methods (Abstract)
        evaluateAt(obj,points)  % evaluation at point list
        fitTo(obj,points) % fit to new data points
    end
    
    % common methods
    methods
       function plot3D(cf,dispfct,dispdim,val_other,ln)         
           np=20;
           stbtemp=cf.SSGrid.StateBounds(:,dispdim);
           totdim=cf.SSGrid.Ndim;
           dim_other=setdiff(1:totdim,dispdim);
           
           gridx=linspace(stbtemp(1,1),stbtemp(2,1),np)';
           gridy=linspace(stbtemp(1,2),stbtemp(2,2),np)';
           
           indmat=grid.StateSpaceGrid.makeCombinations([np,np]);
           
           plist=zeros(np*np,totdim);
           plist(:,dispdim(1))=gridx(indmat(:,2));
           plist(:,dispdim(2))=gridy(indmat(:,1));
           for i=1:length(val_other)
               plist(:,dim_other(i))=ones(np*np,1)*val_other(i);
           end
           
           vlist=evaluateAt(cf,plist);
           vlist=vlist(dispfct,:);
           
           if ln==1
               vlist=exp(vlist);
           elseif ln==2
               vlist=1./(1+exp(vlist));
           end
           figure;
           mesh(gridx,gridy,reshape(vlist,np,np)');
       end        
    end
    
end