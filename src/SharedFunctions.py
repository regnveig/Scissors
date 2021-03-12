from contextlib import contextmanager
from copy import deepcopy as dc
from glob import glob
from multiprocessing import cpu_count, Pool
from pandarallel import pandarallel
from typing import Union
import bz2
import datetime
import functools
import gzip
import io
import json
import logging
import math
import os
import pandas
import pysam
import re
import subprocess
import sys
import tempfile
import time
import warnings

## ------======| LOGGING |======------

def DefaultLogger(
		LogFileName: str,
		Level: int = logging.INFO) -> logging.Logger:
	
	# Format
	Formatter = "%(asctime)-30s%(levelname)-13s%(funcName)-25s%(message)s"
	
	# Compose logger
	Logger = logging.getLogger("default_logger")
	logging.basicConfig(level=Level, format=Formatter)
	
	# Add log file
	Logger.handlers = []
	LogFile = logging.FileHandler(LogFileName)
	LogFile.setLevel(Level)
	LogFile.setFormatter(logging.Formatter(Formatter))
	Logger.addHandler(LogFile)
	
	# Return
	return Logger

## ------======| I/O |======------

def GetContigs(FileBAM: str) -> list: return pysam.AlignmentFile(FileBAM, 'rb').header['SQ']

def SaveJSON(Data: list, FileName: str) -> None: json.dump(Data, open(FileName, 'w'), indent=4, ensure_ascii=False)

def GzipCheck(FileName: str) -> bool: return open(FileName, 'rb').read(2).hex() == "1f8b"

def Bzip2Check(FileName: str) -> bool: return open(FileName, 'rb').read(3).hex() == "425a68"

def OpenAnyway(FileName: str,
		Mode: str,
		Logger: logging.Logger):
	
	try:
		IsGZ = GzipCheck(FileName=FileName)
		IsBZ2 = Bzip2Check(FileName=FileName)
		return gzip.open(FileName, Mode) if IsGZ else (bz2.open(FileName, Mode) if IsBZ2 else open(FileName, Mode))
	except OSError as Err:
		ErrorMessage = f"Can't open the file '{FileName}' ({Err})"
		Logger.error(ErrorMessage)
		raise OSError(ErrorMessage)

def GenerateFileNames(
		Unit: dict,
		Options: dict) -> dict:
	
	Unit['OutputDir'] = os.path.join(Options["PoolDir"], Unit['ID'])
	IRs = os.path.join(Unit['OutputDir'], "IRs")
	FileNames = {
		"OutputDir": Unit['OutputDir'],
		"IRs": IRs,
		"Log": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.pipeline_log.txt"),
		"PrimaryBAM": os.path.join(IRs, f"{Unit['ID']}.primary.bam"),
		"PrimaryStats": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.primary_stats.txt"),
		"DuplessBAM": os.path.join(IRs, f"{Unit['ID']}.dupless.bam"),
		"DuplessMetrics": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.md_metrics.txt"),
		"RecalBAM": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.final.bam"),
		"VCF": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.unfiltered.vcf"),
		"CoverageStats": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.coverage.txt"),
		"VCF": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.unfiltered.vcf"),
		"AnnovarTable": os.path.join(IRs, f"{Unit['ID']}.annovar.tsv"),
		"FilteredXLSX": os.path.join(Unit['OutputDir'], f"{Unit['ID']}.AnnoFit.xlsx")
	}
	return FileNames

## ------======| THREADING |======------

@contextmanager
def Threading(Name: str,
		Logger: logging.Logger,
		Threads: int) -> None:
	
	# Timestamp
	StartTime = time.time()
	
	# Pooling
	pool = Pool(Threads)
	yield pool
	pool.close()
	pool.join()
	del pool
	
	# Timestamp
	Logger.info(f"{Name} finished on {str(Threads)} threads, summary time - %s" % (SecToTime(time.time() - StartTime)))

## ------======| SUBPROCESS |======------

def SimpleSubprocess(
		Name: str,
		Command: str,
		Logger: logging.Logger,
		CheckPipefail: bool = False,
		Env: Union[str, None] = None,
		AllowedCodes: list = []) -> None:
	
	# Timestamp
	StartTime = time.time()
	
	# Compose command
	Command = (f"source {Env}; " if Env is not None else f"") + (f"set -o pipefail; " if CheckPipefail else f"") + Command
	Logger.debug(Command)
	
	# Shell
	Shell = subprocess.Popen(Command, shell=True, executable="/bin/bash", stdout=open(os.devnull, 'w'), stderr=subprocess.PIPE)
	Error = Shell.communicate()[1]
	if Shell.returncode != 0 and Shell.returncode not in AllowedCodes:
		ErrorMessage1 = f"Command '{Name}' has returned non-zero exit code [{str(Shell.returncode)}]."
		ErrorMessage2 = '\n\n' + Error.decode('utf-8') + '\n'
		for line in [ErrorMessage1, Command, ErrorMessage2]: Logger.error(line)
		raise OSError(f"{ErrorMessage1}{ErrorMessage2}")
	if Shell.returncode in AllowedCodes: Logger.warning(f"Command '{Name}' has returned ALLOWED non-zero exit code [{str(Shell.returncode)}].")
	
	# Timestamp
	Logger.info(f"{Name} - %s" % (SecToTime(time.time() - StartTime)))

## ------======| MISC |======------

def SecToTime(Sec: float) -> str: return str(datetime.timedelta(seconds=int(Sec)))

def MultipleTags(Tag: str, List: list, Quoted: bool = True) -> str: return ' '.join([(f"{Tag} \"{str(item)}\"" if Quoted else f"{Tag} {str(item)}") for item in List])

def PrepareGenomeBED(
		Reference: str,
		GenomeBED: str,
		Logger: logging.Logger) -> None:
	
	MODULE_NAME = "PrepareGenomeBED"
	
	# Processing
	SimpleSubprocess(
		Name = f"{MODULE_NAME}.Create",
		Command = "awk 'BEGIN {FS=\"\\t\"}; {print $1 FS \"0\" FS $2}' \"" + Reference + ".fai\" > \"" + GenomeBED + "\"",
		Logger = Logger)

def MakePipe(
		Name: str,
		PipelineConfigFile: str,
		UnitsFile: str) -> tuple:
	
	MODULE_NAME = "MakePipe"
	
	# Options
	CurrentStage = f"{UnitsFile}.{Name}.backup"
	
	# Load data
	PipelineConfig = json.load(open(PipelineConfigFile, 'rt'))
	if os.path.exists(CurrentStage) and os.path.isfile(CurrentStage):
		Protocol = json.load(open(CurrentStage, 'rt'))
		warnings.warn(f"Resume previously interrupted pipeline from backup '{CurrentStage}'")
	else: Protocol = json.load(open(UnitsFile, 'rt'))
	
	# Create backup
	BackupPossible = True
	try:
		SaveJSON(Protocol, CurrentStage)
	except:
		warnings.warn(f"Backup file '{CurrentStage}' cannot be created. If pipeline is interrupted, changes will not be saved")
		BackupPossible = False
	
	# Return
	return (Protocol, PipelineConfig, CurrentStage, BackupPossible)
