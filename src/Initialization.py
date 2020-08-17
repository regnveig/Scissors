from src.SharedFunctions import *

def PrepareReference(Reference, Logger, Env):
	
	MODULE_NAME = "PrepareReference"
	
	# Processing
	SimpleSubprocess(f"{MODULE_NAME}.SamtoolsIndex", f"samtools faidx \"{Reference}\"", Logger)
	SimpleSubprocess(f"{MODULE_NAME}.BWAIndex", f"bwa index \"{Reference}\"", Logger)
	SimpleSubprocess(f"{MODULE_NAME}.GATKIndex",  f"gatk CreateSequenceDictionary -R \"{Reference}\"", Logger, Env=Env) 
