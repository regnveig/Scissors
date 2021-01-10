from contextlib import contextmanager
from glob import glob
from multiprocessing import cpu_count, Pool
from pandarallel import pandarallel
import bz2
import datetime
import functools
import gzip
import json
import logging
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

def DefaultLogger(LogFileName):
	
	# Format
	Formatter = "%(asctime)-30s%(levelname)-13s%(funcName)-25s%(message)s"
	
	# Compose logger
	Logger = logging.getLogger("default_logger")
	logging.basicConfig(level=logging.INFO, format=Formatter)
	
	# Add log file
	Logger.handlers = []
	LogFile = logging.FileHandler(LogFileName)
	LogFile.setLevel(logging.INFO)
	LogFile.setFormatter(logging.Formatter(Formatter))
	Logger.addHandler(LogFile)
	
	# Return
	return Logger


## ------======| SHARED |======------

def GzipCheck(FileName): return open(FileName, 'rb').read(2).hex() == "1f8b"

def Bzip2Check(FileName): return open(FileName, 'rb').read(3).hex() == "425a68"

def SecToTime(Sec): return str(datetime.timedelta(seconds=int(Sec)))

def MultipleTags(Tag, List, Quoted=True): return ' '.join([(f"{Tag} \"{str(item)}\"" if Quoted else f"{Tag} {str(item)}") for item in List])

def GetContigs(FileBAM): return pysam.AlignmentFile(FileBAM, 'rb').header['SQ']

def SaveJSON(Data, FileName): json.dump(Data, open(FileName, 'w'), indent=4, ensure_ascii=False)

def OpenAnyway(FileName, Mode, Logger):
	try:
		IsGZ = GzipCheck(FileName)
		IsBZ2 = Bzip2Check(FileName)
		return gzip.open(FileName, Mode) if IsGZ else (bz2.open(FileName, Mode) if IsBZ2 else open(FileName, Mode))
	except OSError as Err:
		ErrorMessage = f"Can't open the file '{FileName}' ({Err})"
		Logger.error(ErrorMessage)
		raise OSError(ErrorMessage)

@contextmanager
def Threading(Name, Logger, Threads):
	
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

def SimpleSubprocess(Name, Command, Logger, CheckPipefail=False, Env=None, AllowedCodes=[]):
	
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
		Logger.error(ErrorMessage1)
		Logger.error(Command)
		Logger.error(ErrorMessage2)
		raise OSError(f"{ErrorMessage1}{ErrorMessage2}")
	if Shell.returncode in AllowedCodes: Logger.warning(f"Command '{Name}' has returned ALLOWED non-zero exit code [{str(Shell.returncode)}].")
	
	# Timestamp
	Logger.info(f"{Name} - %s" % (SecToTime(time.time() - StartTime))) 
