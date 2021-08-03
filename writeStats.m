% Store all stats in CSV file

% Clear
clc; clear;

% Economies
experdef = '20201117';
suffix = '';
bench_econs = {'base','basegH0','baseC','basesimpleCgH0'};
other_econs = {'theta13','theta14','theta15','theta16','theta17','theta20','theta30',...
                'basegH0','theta14gH0','theta15gH0','theta16gH0',...
                'baseC','theta14C','theta15C','theta16C',...
                'basesimpleCgH0','theta14simpleCgH0','theta15simpleCgH0','theta16simpleCgH0', ...
                'sscpsidn','sscpsiup','sscalphadn','sscalphaup','sscepsilondn','sscepsilonup','sscsigrhoSdn','sscsigrhoSup','piBup','piBdn','basenu','nu10','fullnu2'};
benchecon_index = [1, ones(1,7), 1, 2*ones(1,3), 1, 3*ones(1,3), 1, 4*ones(1,3), ones(1,13)];           
%other_econs = {};

include_all = false;
if include_all
	dircell_sim = struct2cell( dir(['sim_res_',experdef,'_*_',suffix,'.mat']) );
	auto_econs_sim = cellfun(@(x)extractBetween(x,['sim_res_',experdef,'_'],['_',suffix,'.mat']), ...
		dircell_sim(1,:),'UniformOutput',false);
	auto_econs_sim=[auto_econs_sim{:}];
	auto_econs = intersect(auto_econs_sim,auto_econs_irf);
	other_econs = setdiff(auto_econs,{bench_econ});
end
		   
econs = [bench_econs(1),other_econs];

simres_filename = 'simres';
simres_path = ['Results\',simres_filename];
fid=0;

sampleIdentifiers = {'u','e','r','c'};

% Set column headers
colnames={'mean','std','corrY','corrY_1','corrZ','corrZ_1', ...
	'corrVR','corrVR_1','corrGDP','corrGDP_1','AC'};

unitless = {'corrY','corrY_1','corrZ','corrZ_1', ...
	'corrVR','corrVR_1','corrGDP','corrGDP_1','AC'};
isUnitless = ismember(colnames,unitless);

bench_simres = load(['sim_res_',experdef,'_',bench_econs{1},'.mat']);

deffmt = '%s,%0.3f\n';

%% Set parameter formats
paramNames = fieldnames(bench_simres.params);
paramFormats = cell( size( paramNames) );
paramFactors = cell( size( paramNames ) );
paramFormats(:) = {'%0.2f'};
paramFactors(:) = {@(x)x};

match = @(target) cellfun( @(x) contains(x,target), paramNames );
isPct = match('deltaY');
isPct = isPct | strcmp('scaleZ',paramNames);
isPct = isPct | strcmp('Zbar_disc',paramNames);
isPct = isPct | strcmp('piB',paramNames);
isPct = isPct | strcmp('theta',paramNames);
isPct = isPct | strcmp('kappaC',paramNames);
isPct = isPct | strcmp('kappaS',paramNames);
isPct = isPct | strcmp('deltaK',paramNames);
isPct = isPct | strcmp('deltaKbar',paramNames);
isPct = isPct | strcmp('sig_rhoC',paramNames);
isPct = isPct | strcmp('sig_rhoS',paramNames);
isPct = isPct | strcmp('sigZ',paramNames) ; 
isPct = isPct | strcmp('sigY',paramNames) ; 
isPct = isPct | strcmp('xiC',paramNames);
isPct = isPct | strcmp('xiS',paramNames);

isSquare = false;
paramFormats( isPct ) = {'%0.3f'};
paramFactors( isPct & ~isSquare) = {@(x)100*x};
paramFactors( isPct & isSquare) = {@(x)100*sqrt(x)};
paramFactors( ~isPct & isSquare) = {@(x)sqrt(x)};

benchcomp_simres={bench_simres};
for b=2:length(bench_econs)
    benchcomp_simres{b}=load(['sim_res_',experdef,'_',bench_econs{b},'.mat']);
end


cellDict = cell( 1e8, 2 );
ii = 1;
for econCounter = 1:numel(econs)
	econ = econs{econCounter};
	if econCounter > 1
		simres = load(['sim_res_',experdef,'_',econ,'.mat']);
	else
		simres = bench_simres;
		
		% Parameters
		for paramCounter=1:numel(paramNames)
			param = paramNames{paramCounter};
			vecValue = simres.params.(param);
			for paramElementCounter = 1:numel(vecValue)
				paramId = sprintf('param - %s - %d', param,paramElementCounter);
				value = paramFactors{paramCounter}( vecValue( paramElementCounter ) );
				cellDict(ii,:) = writePair(paramId, value, fid, ['%s,',paramFormats{paramCounter},'\n']) ;
				ii = ii + 1;
			end
		end
	end
	
	% Set factors and formats (see function at the bottom of the file)
	[factors, formats] = defineFormats( simres );
	
	% Simulation Results
    for sampleCounter = 1:numel(simres.tabout)
        tab = simres.tabout{sampleCounter};
        thisbench = benchcomp_simres{benchecon_index(econCounter)};
        bench_tab = thisbench.tabout{sampleCounter};
	
        [~,varIdx,varIdx_bench] = intersect( simres.dispnames, ...
            thisbench.dispnames, 'stable' );               
        
        sampleId = sampleIdentifiers{sampleCounter};
        

        for intersectVarCounter = 1:length(varIdx)
            varCounter = varIdx( intersectVarCounter );
            benchVarCounter = varIdx_bench( intersectVarCounter );
            var = simres.dispnames{varCounter};
            fmtStr = formats{varCounter};
            factor = factors{varCounter};
            for statCounter = 1:size(tab,2)
                resId = ['sim - ',econ,' - ',sampleId,' - ',var,' - ',colnames{statCounter}];
                if isUnitless(statCounter)
                    tmpFactor = 1;
                else
                    tmpFactor = factor;
                end
                value = tmpFactor * tab(varCounter,statCounter);
                cellDict(ii,:) = writePair(resId, value, fid, ['%s,',fmtStr,'\n']) ;
                
                compId = ['comp - ',econ,' - ',sampleId,' - ',var,' - ',colnames{statCounter}];
                value = tab(varCounter,statCounter) ./ ...
                    bench_tab(benchVarCounter,statCounter) - 1;
                value = 100 * value;
                cellDict(ii+1,:) = writePair(compId, value, fid, deffmt) ;
                ii = ii + 2;
            end
                        
        end
    end
    
        VHcons = mean( -simres.simseries(:, simres.indexmap.get('VH') ) )^(1-bench_simres.params.gamma);
        cellDict(ii,:) = writePair( ...
			['sim - ',econ,' - u - VHcons - mean'], ...
            VHcons, ...
            fid, deffmt );
        ii = ii+1;
         
        thisbenchsimres=benchcomp_simres{benchecon_index(econCounter)};
        bench_VHcons= mean( -thisbenchsimres.simseries(:, thisbenchsimres.indexmap.get('VH') ) )^(1-thisbenchsimres.params.gamma);
        
         
        cellDict(ii,:) = writePair( ...
            ['comp - ',econ,' - u - VHcons - mean'], ...
            100*(VHcons./bench_VHcons-1), ...
            fid, '%s,%0.2f\n' );
        ii = ii+1;

end



% Process benchmark errors
%errprct = readtable( ['Results\errstats_res_',experdef,'_',bench_econ,'_',suffix,'.xls'] );
errprct = readtable( ['Results\errstats_res_',experdef,'_',bench_econs{1},'fit.xls'] );
labels = errprct.Properties.VariableNames;
for eqnCounter = 1:size(errprct,1)
	for prctCounter = 1:size(errprct,2)
		cellDict(ii,:) = writePair( ...
			sprintf('err - bench - u - eq%02d - %s', eqnCounter, labels{prctCounter} ), ...
			errprct{ eqnCounter, labels{prctCounter} }, ...
			fid, '%s,%0.4f\n' );
		ii = ii + 1;
	end
end

% Wrap up
cellDict = cellDict(1:ii-1,:);
mapDict = containers.Map( cellDict(:,1), cellDict(:,2) );
%fclose(fid);
save([simres_path,'.mat'],'mapDict');



function keyvalcell = writePair( key, val, varargin)
	keyvalcell = { key, val };
	%if nargin > 2
	%	fprintf(fid, fmt, key, val);
	%end
end

% Set variable formats
function [factors, formats, isPct] = defineFormats( simres )
	formats = cell( size( simres.dispnames ) );
	factors = cell( size( simres.dispnames ) );
	formats(:) = {'%0.3f'};
	factors(:) = {1};

	match = @(target) cellfun( @(x) contains(x,target), simres.dispnames );
	isPct = match('log');
	isPct = isPct | match('sh');
    isPct = isPct | match('rate');
	isPct = isPct | match('EER');
	isPct = isPct | match('MRS');
	isPct = isPct | match('lev');
	isPct = isPct | match('FS');
	isPct = isPct | match('FC');
	isPct = isPct | match('DS');
	isPct = isPct | match('DC');
	isPct = isPct | match('recS');
	isPct = isPct | match('recC');
	isPct = isPct | match('wacc');
	isPct = isPct | match('I');
	isPct = isPct | match('assetC');
	isPct = isPct | match('assetS');
	isPct = isPct | match('FinShare');
	isPct = isPct | match('BSsh');
	isPct = isPct | match('KSsh');
    isPct = isPct | match('gYB');
    isPct = isPct | match('liqbenC');
	isPct = isPct | match('AScoeff');
    isPct = isPct | match('ACcoeff');
    isPct = isPct | match('kC');
   % isPct = isPct | match('ACAS');
    

	formats( isPct ) = {'%0.3f'};
	factors( isPct ) = {100};
end
