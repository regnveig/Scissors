from src.SharedFunctions import *
#from src.Regenome import *

# ------======| ANNOVAR |======------

def ANNOVAR(InputVCF, OutputTSV, Databases, DBFolder, AnnovarFolder, GenomeAssembly, Logger, Threads=cpu_count()):
	
	MODULE_NAME = "ANNOVAR"
	
	# Logging
	Logger.info(f"Input VCF: {InputVCF}")
	Logger.info(f"Output TSV: {OutputTSV}")
	Logger.info(f"Genome Assembly: {GenomeAssembly}")
	Logger.info(f"Databases Dir: {DBFolder}")
	Logger.info(f"Databases: {'; '.join([(item['Protocol'] + '[' + item['Operation'] + ']') for item in Databases])}")
	
	# Options
	TableAnnovarPath = os.path.join(AnnovarFolder, "table_annovar.pl")
	Protocol = ','.join([item["Protocol"] for item in Databases])
	Operation = ','.join([item["Operation"] for item in Databases])
	
	with tempfile.TemporaryDirectory() as TempDir:
		TempVCF = os.path.join(TempDir, "temp.vcf")
		AnnotatedTXT = f"{TempVCF}.{GenomeAssembly}_multianno.txt"
		SimpleSubprocess(f"{MODULE_NAME}.TempVCF", f"cp \"{InputVCF}\" \"{TempVCF}\"", Logger)
		SimpleSubprocess(f"{MODULE_NAME}.Annotation", f"perl \"{TableAnnovarPath}\" \"{TempVCF}\" \"{DBFolder}\" --buildver {GenomeAssembly} --protocol {Protocol} --operation {Operation} --remove --vcfinput --thread {Threads}", Logger, AllowedCodes=[25])
		SimpleSubprocess(f"{MODULE_NAME}.CopyTSV", f"cp \"{AnnotatedTXT}\" \"{OutputTSV}\"", Logger)


# ------======| PREPARATION FUNC |======------

def Population(String):
	try:
		return float(String)
	except ValueError:
		return -1.0

def Coords(String):
	try:
		return int(String)
	except ValueError:
		return '.'

def VCF_Data(Block):
	Block["Header"] = [f"VCF.{item}" for item in str(Block["Header"]).split(":")]
	Block["Data"] = str(Block["Data"]).split(":")
	Result = {"Name": Block["Name"]}
	if len(Block["Header"]) != len(Block["Data"]): return Result
	for Num in range(len(Block["Header"])): Result[Block["Header"][Num]] = Block["Data"][Num]
	return Result

def RefMerge(Series):
	Genes = [item for item in Series.to_list() if ((type(item) is str) and (item != "."))]
	Genes = [item.split(';') for item in Genes]
	Genes = list(set([item for sublist in Genes for item in sublist]))
	return ';'.join(sorted(Genes)) if Genes else '.'

def Details(Series):
	lst = [item for item in Series.to_list() if ((type(item) is str) and (item != "."))]
	return '.' if not lst else ';'.join(lst)

def REVEL(Value):
	try:
		return "D" if (float(Value) <= 0.5) else "T"
	except ValueError:
		return "U"

def dbscSNV(Series):
	try:
		return "D" if any([float(item) > 0.6 for item in Series.to_list()]) else "T"
	except ValueError:
		return "."

def MutPred(Str):
	try:
		return "D" if (float(Str) >= 0.9) else "T"
	except ValueError:
		return "U"

def Conservation(Series, Threshold):
	Series = Series.to_list()
	for i in range(len(Series)):
		try:
			Rank = float(Series[i])
			if Rank >= Threshold: Series[i] = "D"
			else: Series[i] = "T"
		except ValueError:
			Series[i] = "U"
	return "/".join([str(Series.count("D")), str(Series.count("D") + Series.count("T"))]) if Series.count("U") < len(Series) else '.'

def SqueezeTable(DataFrame):
	if DataFrame.shape[0] == 1: return DataFrame.iloc[0]
	if DataFrame.shape[0] == 0: return pandas.Series(index=DataFrame.columns.to_list()).fillna('.')
	Squeezed = pandas.Series()
	for col in DataFrame.columns.to_list():
		Squeezed[col] = ';'.join(sorted([str(x) for x in DataFrame[col].to_list() if x != '.']))
		if Squeezed[col] == '': Squeezed[col] = '.'
	return Squeezed

def PredictionsThreshold(Str, Threshold):
	try:
		return int(str(Str).split('/')[0]) >= Threshold
	except ValueError:
		return False

def pLi(Str):
	try:
		return any([float(item) >= 0.9 for item in Str.split(';')])
	except ValueError:
		return False

def Genotype(Str):
	try:
		Lst = [int(item) for item in str(Str).split('/')]
		if len(Lst) != 2: return '.'
		if (Lst[0] == Lst[1]) and (Lst[0] != 0): return 'HOMO'
		else: return Str
	except ValueError:
		return '.'

def Depths(Str):
	try:
		Str = str(Str).split(',')
		return sum([int(item) for item in Str][1:]) >= 4
	except ValueError:
		return False

def FindOMIMCodes(Series):
	Groups = re.findall("\\[MIM:([\\d]+)\\]", str(Series["Disease_description"]))
	Result = {"Name": Series.name}
	for num, item in enumerate(Groups): Result[f"OMIM{str(num)}"] = f"=HYPERLINK(\"https://omim.org/entry/{item}\", \"{item}\")"
	return Result

# ------======| ANNOFIT |======------

def AnnoFit(InputTSV, OutputXLSX, HGMD, AnnovarFolder, AnnoFitConfigFile, Logger, ChunkSize, Threads=cpu_count()):
	
	MODULE_NAME = "AnnoFit"
	
	# Logging
	Logger.info(f"Input TSV: {InputTSV}")
	Logger.info(f"Output XLSX: {OutputXLSX}")
	Logger.info(f"Chunk Size: {str(ChunkSize)}")
	
	GlobalTime = time.time()
	
	# Initialize Pandarallel
	pandarallel.initialize(nb_workers=Threads, verbose=1)
	
	# --- Load Data ---
	
	StartTime = time.time()
	Result = None
	Config = json.load(open(AnnoFitConfigFile, 'rt'))
	HGMDTable = pandas.read_csv(HGMD, sep='\t', dtype=str)
	for Col in ['Chromosome/scaffold position start (bp)', 'Chromosome/scaffold position end (bp)']: HGMDTable[Col] = HGMDTable[Col].parallel_apply(Coords)
	XRefTable = pandas.read_csv(os.path.join(AnnovarFolder, "example/gene_fullxref.txt"), sep='\t', dtype=str).set_index("#Gene_name")
	XRefTable = XRefTable.rename_axis(None, axis=1)
	RegenomeDB = [Regenome(dbname=db["dbName"], filename=db["Filename"], dbtype=db["Type"]) for db in Config["RegenomeAnnotation"]]
	Logger.info(f"Data loaded - %s" % (SecToTime(time.time() - StartTime)))
	
	for ChunkNum, Data in enumerate(pandas.read_csv(InputTSV, sep='\t', dtype=str, chunksize=ChunkSize)):
		
		# --- Regenome ---
		
		StartTime = time.time()
		RDB = Data[["Chr", "Start", "End"]].apply(lambda x: [db.annotation(Data["Chr"], Data["Start"], Data["End"]) for db in RegenomeDB], axis=1)
		print(RDB)
		Logger.info(f"Regenome databases merged - %s" % (SecToTime(time.time() - StartTime)))
		exit()
		
		# --- ANNOVAR Table ---
		
		ChunkTime = time.time()
		StartTime = time.time()
		
		Data.rename(columns=Config["OtherInfo"], inplace=True) # Rename OtherInfo
		Data[["Start", "End"]] = Data[["Start", "End"]].parallel_applymap(Coords) # Prepare coords
		for Col in Config["WipeIntergene"]: Data[Config["WipeIntergene"][Col]] = Data[[Col, Config["WipeIntergene"][Col]]].parallel_apply(lambda x: "." if any([item in Config["IntergeneSynonims"] for item in str(x[Col]).split(';')]) else x[Config["WipeIntergene"][Col]], axis=1)
		Data["AnnoFit.GeneName"] = Data[Config["GeneNames"]].parallel_apply(RefMerge, axis=1) # Gene Names
		Data["AnnoFit.Func"] = Data[Config["Func"]].parallel_apply(RefMerge, axis=1) # Gene Func
		Data["AnnoFit.ExonicFunc"] = Data[Config["ExonicFunc"]].parallel_apply(RefMerge, axis=1) # Gene Exonic Func
		Data["AnnoFit.Details"] = Data[Config["Details"]].parallel_apply(Details, axis=1) # Details
		
		# VCF Data
		VCF_Metadata = pandas.DataFrame(Data[["VCF.FORMAT", "VCF.SAMPLE"]].parallel_apply(lambda x: VCF_Data({"Header": x["VCF.FORMAT"], "Data": x["VCF.SAMPLE"], "Name": x.name}), axis=1).to_list()).set_index("Name")
		del VCF_Metadata.index.name
		Data.drop(columns=["VCF.FORMAT", "VCF.SAMPLE"], inplace=True)
		Data = pandas.concat([Data, VCF_Metadata], axis=1, sort=False)
		del VCF_Metadata
		
		Data["VCF.GT"] = Data["VCF.GT"].parallel_apply(Genotype) # Prepare Genotype
		
		# Symbol Predictions
		for Column in Config["SymbolPred"]: Data[Column["Name"]] = Data[Column["Name"]].parallel_apply(lambda x: "U" if x not in Column["Symbols"] else Column["Symbols"][x])
		Data["REVEL"] = Data["REVEL"].parallel_apply(REVEL)
		Data["MutPred_rankscore"] = Data["MutPred_rankscore"].parallel_apply(MutPred)
		ColNames = [item["Name"] for item in Config["SymbolPred"]] + ["REVEL", "MutPred_rankscore"]
		Data["AnnoFit.ExonPred"] = Data[ColNames].parallel_apply(lambda x: ''.join(x.to_list()), axis=1)
		Data["AnnoFit.ExonPred"] = Data["AnnoFit.ExonPred"].parallel_apply(lambda x: '.' if (x.count('D') + x.count('T') == 0) else '/'.join([str(x.count('D')), str(x.count('D') + x.count('T'))]))
		
		Data["AnnoFit.SplicePred"] = Data[Config["dbscSNV"]].parallel_apply(dbscSNV, axis=1) # dbscSNV
		Data["AnnoFit.Conservation"] = Data[Config["ConservationRS"]].parallel_apply(lambda x: Conservation(x, Threshold=Config['ConservationTH']), axis=1) # Conservation
		
		# Population
		for Col in Config["MedicalPopulationData"]: Data[Col] = Data[Col].parallel_apply(Population)
		Data["AnnoFit.PopFreqMax"] = Data[Config["PopulationData"]].parallel_apply(lambda x: x.apply(Population).max(), axis=1)
		
		Data = Data[Config["ShortVariant"]] # Shorten table
		Logger.info(f"ANNOVAR table is prepared - %s" % (SecToTime(time.time() - StartTime)))
		
		# --- Merge with HGMD ---
		
		StartTime = time.time()
		Data = pandas.merge(Data, HGMDTable, how='left', left_on=["Chr", "Start", "End"], right_on=["Chromosome/scaffold name", "Chromosome/scaffold position start (bp)", "Chromosome/scaffold position end (bp)"])
		Data.rename(columns={"Variant name": "HGMD"}, inplace=True)
		Data["HGMD"].fillna('.', inplace=True)
		Logger.info(f"HGMD merged - %s" % (SecToTime(time.time() - StartTime)))
		
		# --- Merge with XRef ---
		
		StartTime = time.time()
		XRef = Data["AnnoFit.GeneName"].parallel_apply(lambda x: SqueezeTable(XRefTable.loc[[item for item in x.split(';') if item in XRefTable.index],:]))
		Data = pandas.concat([Data, XRef], axis=1, sort=False)
		del XRef
		Logger.info(f"XRef merged - %s" % (SecToTime(time.time() - StartTime)))
		
		# --- Base Filtering ---
		
		StartTime = time.time()
		#DP_filter = Data["VCF.AD"].parallel_apply(Depths)
		#OMIM_filter = Data["Disease_description"].parallel_apply(lambda x: x != '.')
		#HGMD_filter = Data['HGMD'].parallel_apply(lambda x: x != '.')
		#PopMax_filter = Data["AnnoFit.PopFreqMax"] < Config["PopMax_filter"]
		##ExonPred_filter = Data["AnnoFit.ExonPred"].parallel_apply(lambda x: PredictionsThreshold(x, Config["ExonPredThreshold"]))
		#SplicePred_filter = Data["AnnoFit.SplicePred"].parallel_apply(lambda x: x in Config["SplicePred_filter"])
		#IntronPred_filter = Data["regsnp_disease"].parallel_apply(lambda x: x in Config["IntronPred_filter"])
		#Significance_filter = Data["InterVar_automated"].parallel_apply(lambda x: x in Config["InterVar_filter"])
		#CLINVAR_filter = Data["CLNSIG"].parallel_apply(lambda x: any([item in Config["CLINVAR_filter"] for item in str(x).split(',')]))
		##ExonicFunc_filter = Data["AnnoFit.ExonicFunc"].parallel_apply(lambda x: any([item in Config["ExonicFunc_filter"] for item in str(x).split(';')]))
		#ncRNA_filter = Data["AnnoFit.Func"].parallel_apply(lambda x: any([item in Config["ncRNA_filter"] for item in str(x).split(';')]))
		#Splicing_filter = Data["AnnoFit.Func"].parallel_apply(lambda x: any([item in Config["Splicing_filter"] for item in str(x).split(';')]))
		#Data = Data[DP_filter] # & PopMax_filter & (ExonPred_filter | SplicePred_filter | IntronPred_filter | Significance_filter | CLINVAR_filter | ExonicFunc_filter | Splicing_filter | (ncRNA_filter & OMIM_filter))]
		Logger.info(f"Base filtering is ready - %s" % (SecToTime(time.time() - StartTime)))
		
		if Result is None: Result = Data
		else: Result = pandas.concat([Result, Data], axis=0, ignore_index=True)
		
		Logger.info(f"Chunk #{str(ChunkNum + 1)} - %s" % (SecToTime(time.time() - ChunkTime)))

	# Dominance Filtering
	StartTime = time.time()
	#pLi_filter = Result["pLi"].parallel_apply(pLi)
	#OMIM_Dominance_filter = Result["Disease_description"].parallel_apply(lambda x: len(re.findall('[\W]dominant[\W]', str(x).lower())) != 0)
	#Zygocity_filter = Result["VCF.GT"].parallel_apply(lambda x: x == 'HOMO')
	#NoInfo_filter = Result["pLi"].parallel_apply(lambda x: x == '.') & (Result["Disease_description"].parallel_apply(lambda x: (x == '.') | (len(re.findall('([\W]dominant[\W])|([\W]recessive[\W])', str(x).lower())) != 0)))
	Result['Annofit.Compound'] = Result['AnnoFit.GeneName'].parallel_apply(lambda x: ';'.join([str(Result['AnnoFit.GeneName'].apply(lambda y: gene in y.split(';')).value_counts()[True]) for gene in x.split(';')]))
	#Compound_filter = Result['Annofit.Compound'].parallel_apply(lambda x: any([int(item) > 1 for item in x.split(';')]))
	#Result = Result[Compound_filter | pLi_filter | OMIM_Dominance_filter | Zygocity_filter | NoInfo_filter]
	Result = Result.sort_values(by=["Chr", "Start", "End"])
	
	Logger.info(f"Filtering is ready - %s" % (SecToTime(time.time() - StartTime)))
	
	# --- Make Gene list ---
	
	StartTime = time.time()
	Genes = [item.split(';') for item in Result["AnnoFit.GeneName"].to_list()]
	Genes = list(set([item for sublist in Genes for item in sublist]))
	GenesTable = XRefTable.loc[[item for item in Genes if item in XRefTable.index],:]
	GenesTable = GenesTable.reset_index().rename(columns={"index": "#Gene_name"}).sort_values(by=["#Gene_name"])
	GenesTable = GenesTable[Config["ShortGenesTable"]]
	# TODO NoInfo Genes?
	del XRefTable
	Logger.info(f"Genes list is ready - %s" % (SecToTime(time.time() - StartTime)))
	
	# --- Hyperlinks ---
	
	StartTime = time.time()
	Result["avsnp150"] = Result["avsnp150"].parallel_apply(lambda x: '.' if x == '.' else f"=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/snp/{x}\", \"{x}\")")
	Result["UCSC"] = Result[["Chr", "Start", "End"]].parallel_apply(lambda x: f"=HYPERLINK(\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={x['Chr']}%3A{str(x['Start'])}%2D{str(x['End'])}\", \"{x['Chr']}:{str(x['Start'])}\")", axis=1)
	Result["AnnoFit.GeneName"] = Result["AnnoFit.GeneName"].apply(lambda x: '.' if x == '.' else f"=HYPERLINK(\"https://www.genecards.org/Search/Keyword?queryString={'%20OR%20'.join(['%5Baliases%5D(%20' + str(item) + '%20)' for item in x.split(';')])}&keywords={','.join([str(item) for item in x.split(';')])}\", \"{x}\")")
	OMIM_links = pandas.DataFrame(GenesTable[["pLi", "Disease_description"]].parallel_apply(FindOMIMCodes, axis=1).to_list()).set_index("Name").fillna('.')
	del OMIM_links.index.name
	GenesTable = pandas.concat([GenesTable, OMIM_links], axis=1, sort=False)
	Logger.info(f"Hyperlinks are ready - %s" % (SecToTime(time.time() - StartTime)))
	
	# --- Save result ---
	
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

def AnnoPipe(PipelineConfigFile, UnitsFile):
	
	MODULE_NAME = "AnnoPipe"
	
	# Load data
	CurrentStage = f"{UnitsFile}.annotation_stage"
	PipelineConfig = json.load(open(PipelineConfigFile, 'rt'))
	if os.path.exists(CurrentStage) and os.path.isfile(CurrentStage): 
		Units = json.load(open(CurrentStage, 'rt'))
		warnings.warn(f"Resume previously interrupted pipeline from backup '{CurrentStage}'")
	else: Units = json.load(open(UnitsFile, 'rt'))
	
	# Create backup
	BackupPossible = True
	try:
		SaveJSON(Units, CurrentStage)
	except:
		warnings.warn(f"Backup file '{CurrentStage}' cannot be created. If pipeline is interrupted, changes will not be saved")
		BackupPossible = False
	
	# Processing
	for Unit in Units:
		
		# Timestamp
		StartTime = time.time()
		
		# Compose filenames
		IRs = os.path.join(Unit['OutputDir'], "IRs")
		FileNames = {
			"Log": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.pipeline_log.txt"),
			"VCF": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.unfiltered.vcf"),
			"AnnovarTable": os.path.join(IRs, f"{Unit['ID']}.annovar.tsv"),
			"FilteredXLSX": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.AnnoFit.xlsx")
			}
		
		# Configure logging
		Logger = DefaultLogger(FileNames["Log"])
		
		if Unit["Stage"] == 0:
			
			# ANNOVAR
			ANNOVAR(FileNames["VCF"], FileNames["AnnovarTable"], PipelineConfig["AnnovarDatabases"], PipelineConfig["AnnovarDBFolder"], PipelineConfig["AnnovarFolder"], PipelineConfig["GenomeAssembly"], Logger, Threads=PipelineConfig["Threads"])
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
		
		if Unit["Stage"] == 1:
			
			# AnnoFit
			AnnoFit(FileNames["AnnovarTable"], FileNames["FilteredXLSX"], PipelineConfig["HGMDPath"], PipelineConfig["AnnovarFolder"], PipelineConfig["AnnoFitConfig"], Logger, PipelineConfig["AnnofitChunkSize"], Threads=PipelineConfig["Threads"])
			
			# Timestamp
			Logger.info(f"Unit {Unit['ID']} successfully annotated, summary time - %s" % (SecToTime(time.time() - StartTime)))
			
			Unit["Stage"] += 1
			if BackupPossible: SaveJSON(Units, CurrentStage)
	
	os.remove(CurrentStage)
