function status = build(texfile,dbfile,compiler)

optionalDefaults = {'mean','u','bench','2'};
paramOptionalDefaults = {'1'};

%pattern{1} = '\\simres\{(?<var>.*?)\}(?<stat>\[.*?\])?(?<sample>\[.*?\])?(?<econ>\[.*?\])?(?<round>\[.*?\])?';
%pattern{1} = '\\simres\{(?<var>.*?)\}(\[(?<stat>.*?)\])';
pattern = '\\(?<set>[a-z]+?)res\{(?<var>.*?)\}(?<optional>(\[.*?\]){0,4})';
paramPattern = '\\param\{(?<param>.*?)\}(?<optional>(\[.*?\])?)';

load(dbfile,'mapDict');
paper = fileread([texfile,'.tex']);

paramStr = 'param';
replace = '${replaceVal($2,$3,$1,$0,mapDict,optionalDefaults)}';
paramReplace = '${replaceVal($1,$2,paramStr,$0,mapDict,paramOptionalDefaults)}';
newpaper = regexprep( paper, pattern, replace );
newpaper = regexprep( newpaper, paramPattern, paramReplace );
%newpaper = regexprep( newpaper, '\\DTLloaddbtex', '%\\DTLloaddbtex');

fid = fopen([texfile,'_filled.tex'],'w');
fprintf(fid,'%s',newpaper);
fclose(fid);

if ~isempty(compiler)
	[tok,arr] = regexp(texfile,'([\\/])','tokens','split');
	if length(arr)>1
		folder = strjoin(arr(1:end-1),tok{1}{1});
	else
		folder = cd;
	end
	file = arr{end};
	
	currdir = cd(folder);
	compileStr = [compiler,'.exe -synctex=1 -interaction=nonstopmode ', ...
		file,'_filled.tex"'];
	bibStr = ['bibtex.exe "',file,'_filled"'];	
	
	try
		out(1) = system(compileStr);
		out(2) = system(bibStr);
		out(3) = system(compileStr);
		out(4) = system(compileStr);
		status = ~all(out==0);
	catch e
		cd(currdir);
		status = -1;
		rethrow(e);
	end
	copyfile([file,'_filled.pdf'],[file,'.pdf']);
	%delete([file,'_filled.*']);
	cd(currdir);
end

end
