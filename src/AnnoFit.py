from src.SharedFunctions import *

# ------======| ANNOVAR |======------

def ANNOVAR(
		InputVCF: str,
		OutputTSV: str,
		DBFolder: str,
		AnnovarFolder: str,
		GenomeAssembly: str,
		Logger: logging.Logger,
		Databases: list = [],
		GFF3List: list = [],
		Threads: int = cpu_count()) -> None:
	
	MODULE_NAME = "ANNOVAR"
	
	assert bool(Databases) != bool(GFF3List), f"Either Databases or GFF3 files must be defined"
	
	# Logging
	for line in [f"Input VCF: {InputVCF}", f"Output TSV: {OutputTSV}", f"Genome Assembly: {GenomeAssembly}", f"Databases Dir: {DBFolder}"] + ([] if not Databases else [f"Databases: {'; '.join([(item['Protocol'] + '[' + item['Operation'] + ']') for item in Databases])}"]) + ([] if not GFF3List else [f"Databases: GFF3, {len(GFF3List)} items [r]"]): Logger.info(line)
	
	with tempfile.TemporaryDirectory() as TempDir:
		
		# Options
		TableAnnovarPath = os.path.join(AnnovarFolder, "table_annovar.pl")
		Protocol = ','.join([item["Protocol"] for item in Databases] + ["gff3" for item in GFF3List])
		Operation = ','.join([item["Operation"] for item in Databases] + ["r" for item in GFF3List])
		GFFs = "--gff3dbfile " + ','.join(GFF3List)
		TempVCF = os.path.join(TempDir, "temp.vcf")
		AnnotatedTXT = f"{TempVCF}.{GenomeAssembly}_multianno.txt"
		
		# Processing
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.TempVCF",
			Command = f"cp \"{InputVCF}\" \"{TempVCF}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.Annotation",
			Command = f"perl \"{TableAnnovarPath}\" \"{TempVCF}\" \"{DBFolder}\" --buildver {GenomeAssembly} --protocol {Protocol} --operation {Operation} {GFFs} --remove --vcfinput --thread {Threads}",
			Logger = Logger,
			AllowedCodes = [25])
		SimpleSubprocess(
			Name = f"{MODULE_NAME}.CopyTSV",
			Command = f"cp \"{AnnotatedTXT}\" \"{OutputTSV}\"",
			Logger = Logger)

# ------======| CUSTOM REGION-BASED ANNOTATIONS |======------

def Tsv2Gff3(
		dbName: str,
		InputTSV: str,
		ChromCol: str,
		StartCol: str,
		EndCol: str,
		OutputGFF3: str,
		Reference: str,
		Logger: logging.Logger,
		Threads: int) -> None:
	
	MODULE_NAME = "Tsv2Gff3"
	
	pandarallel.initialize(nb_workers=Threads, verbose=1)
	
	# Logging
	for line in [f"Name: {dbName}", f"Input TSV db: {InputTSV}", f"Output GFF3: {OutputGFF3}"]: Logger.info(line)
	
	# Options
	AnchorCols = [ChromCol, StartCol, EndCol]
	DataOrder = [ChromCol, "sample", "type", StartCol, EndCol, "score", "strand", "phase", "attributes"]
	
	# Prepare faidx
	Faidx = pandas.read_csv(Reference + ".fai", sep='\t', header=None).assign(Tag="##sequence-region", Start=1)[["Tag", 0, "Start", 1]]
	Chroms = {value: index for index, value in enumerate(Faidx[0].to_list())}
	
	# Load data
	Data = pandas.read_csv(InputTSV, sep='\t')
	AttributeCols = [item for item in Data.columns.to_list() if item not in AnchorCols]
	
	# Filter & sort intervals by reference
	Filtered = [item for item in list(set(Data[ChromCol].to_list())) if item not in Chroms.keys()]
	if Filtered: Logger.warning(f"Contigs will be removed from database \"{dbName}\": {', '.join(sorted(Filtered))}")
	Data = Data[Data[ChromCol].parallel_apply(lambda x: x in Chroms.keys())]
	Data["Rank"] = Data[ChromCol].map(Chroms)
	Data.sort_values(["Rank", StartCol], inplace=True)
	if Data.shape[0] == 0:
		ErrorMessage = f"Database \"{dbName}\" and reference \"{Reference}\" have no matching columns"
		Logger.error(ErrorMessage)
		raise RuntimeError(ErrorMessage)
	
	# Processing
	Attributes = Data[AttributeCols].parallel_apply(lambda x: "ID=" + base64.b16encode(json.dumps(x.to_dict()).encode('utf-8')).decode('utf-8'), axis=1)
	Data.drop(columns=AttributeCols, inplace=True)
	Data[StartCol] = Data[StartCol].parallel_apply(lambda x: 1 if x == 0 else x)
	Data = Data.assign(sample=dbName, type="region", attributes=Attributes,	score=".", strand=".", phase=".")[DataOrder]
	
	# Save
	with open(OutputGFF3, 'wt') as O:
		O.write(
			"##gff-version 3\n" +
			Faidx.to_csv(sep=' ', index=False, header=False) + 
			Data.to_csv(sep='\t', index=False, header=False))
	
	# Return expected cols
	return [f"{str(dbName)}.{str(item)}" for item in AttributeCols] if AttributeCols else [ str(dbName) ]

def CureBase(
		InputVCF: str,
		OutputTSV: str,
		Databases: list,
		AnnovarFolder: str,
		GenomeAssembly: str,
		Reference: str,
		Logger: logging.Logger,
		Threads: int) -> None:
	
	MODULE_NAME = "CureBase"
	
	# Initialize Pandarallel
	pandarallel.initialize(nb_workers=Threads, verbose=1)
	
	with tempfile.TemporaryDirectory() as TempDir:
		
		Gff3List, ExpectedCols = [], []
		
		for index, DB in enumerate(Databases):
			Gff3File = os.path.join(TempDir, f"database_{str(index)}.gff3")
			ExpectedCols += Tsv2Gff3(
				dbName = DB["Name"],
				InputTSV = DB["FileName"],
				ChromCol = DB["ChromColumn"],
				StartCol = DB["StartColumn"],
				EndCol = DB["EndColumn"],
				OutputGFF3 = Gff3File,
				Reference = Reference,
				Logger = Logger,
				Threads = Threads)
			Gff3List += [ f"database_{str(index)}.gff3" ]
		
		TempTSV = os.path.join(TempDir, f"temp.tsv")
		ANNOVAR(
			InputVCF = InputVCF,
			OutputTSV = TempTSV,
			GFF3List = Gff3List,
			DBFolder = TempDir,
			AnnovarFolder = AnnovarFolder,
			GenomeAssembly = GenomeAssembly,
			Logger = Logger,
			Threads = Threads)
		
		SNPdata = ['Chr', 'Start', 'End', 'Ref', 'Alt']
		Data = pandas.read_csv(TempTSV, sep='\t', dtype=str)
		NewColumns = {f"gff3{'' if index == 0 else str(index + 1)}": item["Name"] for index, item in enumerate(Databases)}
		Data = Data[SNPdata + list(NewColumns.keys())]
		Data[list(NewColumns.keys())] = Data[list(NewColumns.keys())].parallel_applymap(lambda x: '.' if x == '.' else [json.loads(base64.b16decode(item.encode('utf-8')).decode('utf-8')) for item in x.split("=")[1].split(",")])
		for Col in list(NewColumns.keys()):
			NewCols = Data[Col][Data[Col] != '.']
			if NewCols.size > 0:
				NewCols = NewCols.parallel_apply(lambda LD: pandas.Series({str(NewColumns[Col]): "yes"} if not LD[0] else {f"{str(NewColumns[Col])}.{str(k)}": [dic[k] for dic in LD] for k in LD[0]}))
				Data = pandas.concat([Data, NewCols], axis=1)
			else: Logger.warning(f"Database has no intersections with variants: {str(NewColumns[Col])}")
		Data = Data.drop(columns=list(NewColumns.keys()))
		MissingCols = [item for item in ExpectedCols if item not in Data.columns.to_list()]
		Data[MissingCols] = float("nan")
		InformationCols = [item for item in Data.columns.to_list() if item not in SNPdata]
		Data[InformationCols] = Data[InformationCols].parallel_applymap(lambda x: '.' if x != x else '; '.join([str(item) for item in sorted(list(set(x)))]))
		Data = Data[SNPdata + ExpectedCols]
		Data.to_csv(OutputTSV, sep='\t', index=False)

# ------======| ANNOFIT |======------

def AnnoFit(
		InputTSV: str,
		OutputXLSX: str,
		HGMD: str,
		AnnovarFolder: str,
		AnnoFitConfigFile: str,
		Logger: logging.Logger,
		ChunkSize: int,
		Threads: int = cpu_count()) -> None:
	
	MODULE_NAME = "AnnoFit"
	
	# Logging
	for line in [f"Input TSV: {InputTSV}", f"Output XLSX: {OutputXLSX}", f"Chunk Size: {str(ChunkSize)}"]: Logger.info(line)
	
	# Initialize Pandarallel
	pandarallel.initialize(nb_workers=Threads, verbose=1)
	
	# Global func
	def FormatInt(String: str) -> Union[int, None]:
		try:
			return int(String)
		except ValueError:
			return None
	def FormatFloat(String: str) -> Union[float, None]:
		try:
			return float(String)
		except ValueError:
			return None
	def SqueezeTable(DataFrame: pandas.DataFrame) -> pandas.Series:
		if DataFrame.shape[0] == 1: return DataFrame.iloc[0]
		if DataFrame.shape[0] == 0: return pandas.Series(index=DataFrame.columns.to_list()).fillna('.')
		Squeezed = pandas.Series()
		for col in DataFrame.columns.to_list():
			Squeezed[col] = ';'.join(sorted([str(x) for x in DataFrame[col].to_list() if x != '.']))
			if Squeezed[col] == '': Squeezed[col] = '.'
		return Squeezed
	
	# Format func
	def FormatCoordinates(String: str) -> Union[int, str]:
		Result = FormatInt(String)
		return ('.' if Result is None else Result)
	def FormatPopulationFreq(String: str) -> float:
		Result = FormatFloat(String)
		return (-1.0 if Result is None else Result)
	def FormatVcfMetadata(Block: pandas.Series) -> dict:
		Block = {"Header": Block["VCF.FORMAT"], "Data": Block["VCF.SAMPLE"], "Name": Block.name}
		Block["Header"] = [f"VCF.{item}" for item in str(Block["Header"]).split(":")]
		Block["Data"] = str(Block["Data"]).split(":")
		Result = {"Name": Block["Name"]}
		if len(Block["Header"]) != len(Block["Data"]): return Result
		for Num in range(len(Block["Header"])): Result[Block["Header"][Num]] = Block["Data"][Num]
		return Result
	def FormatGenesOrFunction(Series: pandas.Series) -> str:
		Genes = [item for item in Series.to_list() if ((type(item) is str) and (item != "."))]
		Genes = [item.split(';') for item in Genes]
		Genes = list(set([item for sublist in Genes for item in sublist]))
		return (';'.join(sorted(Genes)) if Genes else '.')
	def FormatRevel(Value: str) -> str:
		Value = FormatFloat(Value)
		return ("U" if Value is None else ("D" if (Value <= 0.5) else "T"))
	def FormatDbscSNV(Series: pandas.Series) -> str:
		Result = [FormatFloat(item) for item in Series.to_list()]
		return ("." if any([item is None for item in Result]) else ("D" if any([item > 0.6 for item in Result]) else "T"))
	def FormatMutPred(Str: str) -> str:
		Result = FormatFloat(Str)
		return ("U" if Result is None else ("D" if (Result >= 0.9) else "T"))
	def FormatOmimCodes(Series: pandas.Series) -> dict:
		Groups = re.findall("\\[MIM:([\\d]+)\\]", str(Series["Disease_description"]))
		Result = {"Name": Series.name}
		for num, item in enumerate(Groups): Result[f"OMIM-{num:02d}"] = f"=HYPERLINK(\"https://omim.org/entry/{item}\", \"{item}\")"
		return Result
	def FormatGenotype(Str: str) -> str:
		Lst = [FormatInt(item) for item in Str.split('/')]
		if (len(Lst) != 2) or any([item is None for item in Lst]): return '.'
		return ('HOMO' if ((Lst[0] == Lst[1]) and (Lst[0] != 0)) else Str)
	def FormatConservation(Series: pandas.Series) -> str:
		Series = Series.to_list()
		for i in range(len(Series)):
			Rank = FormatFloat(Series[i])
			Series[i] = "U" if Rank is None else ("D" if Rank >= 0.7 else "T")
		return "/".join([str(Series.count("D")), str(Series.count("D") + Series.count("T"))]) if Series.count("U") < len(Series) else '.'
	def FormatDetails(Series: pandas.Series) -> str:
		lst = [item for item in Series.to_list() if ((type(item) is str) and (item != "."))]
		return '.' if not lst else ';'.join(lst)
	def FormatExonPrediction(Series: pandas.Series) -> str:
		Result = ''.join(Series.to_list())
		return ('.' if (Result.count('D') + Result.count('T') == 0) else '/'.join([str(Result.count('D')), str(Result.count('D') + Result.count('T'))]))
	FormatdbSNP = lambda x: '.' if x == '.' else f"=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/snp/{x}\", \"{x}\")"
	FormatUCSC = lambda x: f"=HYPERLINK(\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={x['Chr']}%3A{str(x['Start'])}%2D{str(x['End'])}\", \"{x['Chr']}:{str(x['Start'])}\")"
	FormatGenomeBrowser = lambda x: '.' if x == '.' else f"=HYPERLINK(\"https://www.genecards.org/Search/Keyword?queryString={'%20OR%20'.join(['%5Baliases%5D(%20' + str(item) + '%20)' for item in x.split(';')])}&keywords={','.join([str(item) for item in x.split(';')])}\", \"{x}\")"
	
	# Filter func
	def FilterpLi(Str: str) -> bool:
		Result = [FormatFloat(item) for item in Str.split(';')]
		return (False if any([item is None for item in Result]) else any([item >= 0.9 for item in Result]))
	def FilterDepth(Str: str) -> bool:
		Result = [FormatInt(item) for item in Str.split(',')]
		return (False if any([item is None for item in Result]) else any([item >= 4 for item in Result]))
	def FilterExonPrediction(Str: str) -> bool:
		Result = FormatInt(Str.split('/')[0])
		return (False if Result is None else (Result >= 3))
	FilterOmimDominance = lambda x: len(re.findall('[\W]dominant[\W]', str(x).lower())) != 0
	FilterNoInfo = lambda x: (x["pLi"] == '.') and (x["Disease_description"] == '.') or (len(re.findall('([\W]dominant[\W])|([\W]recessive[\W])', str(["Disease_description"]).lower())) == 0)
	
	# Load Data
	GlobalTime = time.time()
	StartTime = time.time()
	Result = None
	Config = json.load(open(AnnoFitConfigFile, 'rt'))
	HGMDTable = pandas.read_csv(HGMD, sep='\t', dtype=str)
	for Col in ['Chromosome/scaffold position start (bp)', 'Chromosome/scaffold position end (bp)']: HGMDTable[Col] = HGMDTable[Col].parallel_apply(FormatCoordinates)
	XRefTable = pandas.read_csv(os.path.join(AnnovarFolder, "example/gene_fullxref.txt"), sep='\t', dtype=str).set_index("#Gene_name").rename_axis(None, axis=1)
	Logger.info(f"Data loaded - %s" % (SecToTime(time.time() - StartTime)))
	
	for ChunkNum, Data in enumerate(pandas.read_csv(InputTSV, sep='\t', dtype=str, chunksize=ChunkSize)):
		
		# ANNOVAR Table
		ChunkTime = time.time()
		StartTime = time.time()
		
		# Basic info format
		Data.rename(columns=Config["OtherInfo"], inplace=True) # Rename OtherInfo
		Data[["Start", "End"]] = Data[["Start", "End"]].parallel_applymap(FormatCoordinates) # Prepare coords
		for Col in Config["WipeIntergene"]: Data[Config["WipeIntergene"][Col]] = Data[[Col, Config["WipeIntergene"][Col]]].parallel_apply(lambda x: "." if any([item in Config["IntergeneSynonims"] for item in str(x[Col]).split(';')]) else x[Config["WipeIntergene"][Col]], axis=1)
		Data["AnnoFit.GeneName"] = Data[Config["GeneNames"]].parallel_apply(FormatGenesOrFunction, axis=1) # Gene Names
		Data["AnnoFit.Func"] = Data[Config["Func"]].parallel_apply(FormatGenesOrFunction, axis=1) # Gene Func
		Data["AnnoFit.ExonicFunc"] = Data[Config["ExonicFunc"]].parallel_apply(FormatGenesOrFunction, axis=1) # Gene Exonic Func
		Data["AnnoFit.Details"] = Data[Config["Details"]].parallel_apply(FormatDetails, axis=1) # Details
		
		# VCF Data
		VCF_Metadata = pandas.DataFrame(Data[["VCF.FORMAT", "VCF.SAMPLE"]].parallel_apply(FormatVcfMetadata, axis=1).to_list()).set_index("Name").rename_axis(None, axis=1)
		Data.drop(columns=["VCF.FORMAT", "VCF.SAMPLE"], inplace=True)
		Data = pandas.concat([Data, VCF_Metadata], axis=1, sort=False)
		del VCF_Metadata
		Data["VCF.GT"] = Data["VCF.GT"].parallel_apply(FormatGenotype) # Prepare Genotype
		
		# Symbol Predictions
		for Column in Config["SymbolPred"]: Data[Column["Name"]] = Data[Column["Name"]].parallel_apply(lambda x: "U" if x not in Column["Symbols"] else Column["Symbols"][x])
		Data["REVEL"] = Data["REVEL"].parallel_apply(FormatRevel)
		Data["MutPred_rankscore"] = Data["MutPred_rankscore"].parallel_apply(FormatMutPred)
		ColNames = [item["Name"] for item in Config["SymbolPred"]] + ["REVEL", "MutPred_rankscore"]
		Data["AnnoFit.ExonPred"] = Data[ColNames].parallel_apply(FormatExonPrediction, axis=1)
		Data["AnnoFit.SplicePred"] = Data[Config["dbscSNV"]].parallel_apply(FormatDbscSNV, axis=1) # dbscSNV
		Data["AnnoFit.Conservation"] = Data[Config["ConservationRS"]].parallel_apply(FormatConservation, axis=1) # Conservation
		
		# Population
		for Col in Config["MedicalPopulationData"]: Data[Col] = Data[Col].parallel_apply(FormatPopulationFreq)
		Data["AnnoFit.PopFreqMax"] = Data[Config["PopulationData"]].parallel_apply(lambda x: x.apply(FormatPopulationFreq).max(), axis=1)
		
		# Shorten table
		Data = Data[Config["ShortVariant"]]
		Logger.info(f"ANNOVAR table is prepared - %s" % (SecToTime(time.time() - StartTime)))
		
		# Merge with HGMD
		StartTime = time.time()
		Data = pandas.merge(Data, HGMDTable, how='left', left_on=["Chr", "Start", "End"], right_on=["Chromosome/scaffold name", "Chromosome/scaffold position start (bp)", "Chromosome/scaffold position end (bp)"])
		Data.rename(columns={"Variant name": "HGMD"}, inplace=True)
		Data["HGMD"].fillna('.', inplace=True)
		Logger.info(f"HGMD merged - %s" % (SecToTime(time.time() - StartTime)))
		
		# Merge with XRef
		StartTime = time.time()
		XRef = Data["AnnoFit.GeneName"].parallel_apply(lambda x: SqueezeTable(XRefTable.loc[[item for item in x.split(';') if item in XRefTable.index],:]))
		Data = pandas.concat([Data, XRef], axis=1, sort=False)
		del XRef
		Logger.info(f"XRef merged - %s" % (SecToTime(time.time() - StartTime)))
		
		# Base Filtering
		StartTime = time.time()
		Filters = {}
		Filters["DP"] = Data["VCF.AD"].parallel_apply(FilterDepth)
		Filters["OMIM"] = Data["Disease_description"].parallel_apply(lambda x: x != '.')
		Filters["HGMD"] = Data['HGMD'].parallel_apply(lambda x: x != '.')
		Filters["PopMax"] = Data["AnnoFit.PopFreqMax"] < Config["PopMax_filter"]
		Filters["ExonPred"] = Data["AnnoFit.ExonPred"].parallel_apply(FilterExonPrediction)
		Filters["SplicePred"] = Data["AnnoFit.SplicePred"].parallel_apply(lambda x: x in Config["SplicePred_filter"])
		Filters["IntronPred"] = Data["regsnp_disease"].parallel_apply(lambda x: x in Config["IntronPred_filter"])
		Filters["Significance"] = Data["InterVar_automated"].parallel_apply(lambda x: x in Config["InterVar_filter"])
		Filters["CLINVAR"] = Data["CLNSIG"].parallel_apply(lambda x: any([item in Config["CLINVAR_filter"] for item in str(x).split(',')]))
		Filters["ExonicFunc"] = Data["AnnoFit.ExonicFunc"].parallel_apply(lambda x: any([item in Config["ExonicFunc_filter"] for item in str(x).split(';')]))
		Filters["ncRNA"] = Data["AnnoFit.Func"].parallel_apply(lambda x: any([item in Config["ncRNA_filter"] for item in str(x).split(';')]))
		Filters["Splicing"] = Data["AnnoFit.Func"].parallel_apply(lambda x: any([item in Config["Splicing_filter"] for item in str(x).split(';')]))
		#Data = Data[Filters["DP"]] # & PopMax_filter & (ExonPred_filter | SplicePred_filter | IntronPred_filter | Significance_filter | CLINVAR_filter | ExonicFunc_filter | Splicing_filter | (ncRNA_filter & OMIM_filter))]
		Logger.info(f"Base filtering is ready - %s" % (SecToTime(time.time() - StartTime)))
		
		#Concat chunks
		Result = Data if Result is None else pandas.concat([Result, Data], axis=0, ignore_index=True)
		Logger.info(f"Chunk #{str(ChunkNum + 1)} - %s" % (SecToTime(time.time() - ChunkTime)))
	
	# Compound
	Result['Annofit.Compound'] = Result['AnnoFit.GeneName'].parallel_apply(lambda x: ';'.join([str(Result['AnnoFit.GeneName'].apply(lambda y: gene in y.split(';')).value_counts()[True]) for gene in x.split(';')]))
	
	# Dominance Filtering
	StartTime = time.time()
	Filters["pLi"] = Result["pLi"].parallel_apply(FilterpLi)
	Filters["OMIM_Dominance"] = Result["Disease_description"].parallel_apply(FilterOmimDominance)
	Filters["Zygocity"] = Result["VCF.GT"].parallel_apply(lambda x: x == 'HOMO')
	Filters["NoInfo"] = Result[["pLi", "Disease_description"]].parallel_apply(FilterNoInfo, axis=1)
	Filters["Compound_filter"] = Result['Annofit.Compound'].parallel_apply(lambda x: any([int(item) > 1 for item in x.split(';')]))
	#Result = Result[Compound_filter | pLi_filter | OMIM_Dominance_filter | Zygocity_filter | NoInfo_filter]
	Result = Result.sort_values(by=["Chr", "Start", "End"])
	Logger.info(f"Filtering is ready - %s" % (SecToTime(time.time() - StartTime)))
	
	# Make genes list
	StartTime = time.time()
	Genes = [item.split(';') for item in Result["AnnoFit.GeneName"].to_list()]
	Genes = list(set([item for sublist in Genes for item in sublist]))
	GenesTable = XRefTable.loc[[item for item in Genes if item in XRefTable.index],:].reset_index().rename(columns={"index": "#Gene_name"}).sort_values(by=["#Gene_name"])[Config["ShortGenesTable"]]
	# TODO NoInfo Genes?
	del XRefTable
	Logger.info(f"Genes list is ready - %s" % (SecToTime(time.time() - StartTime)))
	
	# Hyperlinks
	StartTime = time.time()
	Result["avsnp150"] = Result["avsnp150"].parallel_apply(FormatdbSNP)
	Result["UCSC"] = Result[["Chr", "Start", "End"]].parallel_apply(FormatUCSC, axis=1)
	Result["AnnoFit.GeneName"] = Result["AnnoFit.GeneName"].apply(FormatGenomeBrowser)
	OMIM_links = pandas.DataFrame(GenesTable[["pLi", "Disease_description"]].parallel_apply(FormatOmimCodes, axis=1).to_list()).set_index("Name").fillna('.').rename_axis(None, axis=1)
	GenesTable = pandas.concat([GenesTable, OMIM_links], axis=1, sort=False)
	Logger.info(f"Hyperlinks are ready - %s" % (SecToTime(time.time() - StartTime)))
	
	# Save result
	StartTime = time.time()
	Result = Result[Config["FinalVariant"]]
	Result.insert(0, 'Comment', '')
	with pandas.ExcelWriter(OutputXLSX) as Writer:
		Result.to_excel(Writer, "Variants", index=False)
		GenesTable.to_excel(Writer, "Genes", index=False)
		Writer.save()
	Logger.info(f"Files saved - %s" % (SecToTime(time.time() - StartTime)))
	Logger.info(f"{MODULE_NAME} finish - %s" % (SecToTime(time.time() - GlobalTime)))

# ------======| ANNOTATION PIPELINE |======------

def AnnoPipe(
		PipelineConfigFile: str,
		UnitsFile: str) -> None:
	
	MODULE_NAME = "AnnoPipe"
	Protocol, PipelineConfig, CurrentStage, BackupPossible = MakePipe(MODULE_NAME, PipelineConfigFile, UnitsFile) # MakePipe
	
	# Processing
	for Unit in Protocol["Units"]:
		StartTime = time.time()
		FileNames = GenerateFileNames(Unit, Protocol["Options"]) # Compose filenames
		Logger = DefaultLogger(FileNames["Log"]) # Configure logging
		
		if Unit["Stage"] == 0:
			
			ANNOVAR(
				InputVCF = FileNames["VCF"],
				OutputTSV = FileNames["AnnovarTable"],
				Databases = PipelineConfig["AnnovarDatabases"],
				DBFolder = PipelineConfig["AnnovarDBFolder"],
				AnnovarFolder = PipelineConfig["AnnovarFolder"],
				GenomeAssembly = PipelineConfig["GenomeAssembly"],
				Logger = Logger,
				Threads = PipelineConfig["Threads"])
			
			if PipelineConfig["GFF3"]:
				CureBase(
					InputVCF = FileNames["VCF"],
					OutputTSV = FileNames["Gff3Table"],
					Databases = PipelineConfig["GFF3"],
					AnnovarFolder = PipelineConfig["AnnovarFolder"],
					Reference = PipelineConfig["Reference"],
					GenomeAssembly = PipelineConfig["GenomeAssembly"],
					Logger = Logger,
					Threads = PipelineConfig["Threads"])
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
		
		if Unit["Stage"] == 2:
			AnnoFit(
				InputTSV = FileNames["AnnovarTable"],
				OutputXLSX = FileNames["FilteredXLSX"],
				HGMD = PipelineConfig["HGMDPath"],
				AnnovarFolder = PipelineConfig["AnnovarFolder"],
				AnnoFitConfigFile = PipelineConfig["AnnoFitConfig"],
				Logger = Logger,
				ChunkSize = PipelineConfig["AnnofitChunkSize"],
				Threads = PipelineConfig["Threads"])
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Protocol, CurrentStage)
			Logger.info(f"Unit {Unit['ID']} successfully annotated, summary time - %s" % (SecToTime(time.time() - StartTime)))
	
	os.remove(CurrentStage)
