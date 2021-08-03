classdef HelperCollection
   
    properties
       dummy         
    end
    
    methods
        function obj=HelperCollection()
            obj.dummy=1;
            disp('No need to create an instance of this class.');
        end
    end
    
    methods (Static)
        
        function bigstr=combineStructs(strArray)
            % strArray must be cell array of structs
            
            bigstr=struct;
            Ninstr=length(strArray);
            for i=1:Ninstr
                strtmp=strArray{i};
                fntmp=fieldnames(strtmp);
                for j=1:length(fntmp)
                    thisfield=fntmp{j};
                    bigstr.(thisfield)=strtmp.(thisfield);
                end
            end
        end
        
        function [indlist,effnames]=makeListFromNames(hashmap,names)            
            n=length(names);
            indlist=zeros(n,1);
            for i=1:n
                tmp=hashmap.get(names{i});
                if ~isempty(tmp)
                    indlist(i)=tmp;
                else
                    indlist(i)=0;
                end
            end
            effnames=names(indlist~=0);
            indlist=indlist(indlist~=0);            
        end
        
        function dum=tableExport(filename, varnames, data)
            
            dum=[];
            
            ncol=length(varnames);
            fid=fopen(filename,'w');
            for c=1:ncol-1
                fprintf(fid,'%s,',varnames{c});
            end
            fprintf(fid,'%s\n',varnames{ncol});
            fclose(fid);
            
            dlmwrite(filename,data,'-append');
            
        end
        
        function dum=scatterPoints2D(X,Y)
            dum=[];
            
            bY=max(Y);
            scale=100/bY;
            figure; hold on;
            scatter(X(:,1),X(:,2),Y*scale);            
        end
        
    end
    
end