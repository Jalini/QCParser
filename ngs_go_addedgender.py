#!/usr/bin/env python
import argparse
import time
import datetime
import sys
import dragen
import os
import picard_util
import qdx_util
import re
import shutil
import subprocess


#pipeline to run process from bcl to database
#jje, ce 01242019
#r&d@qdx
#three parts: bcl2fastq, variant calling, and qc report
iscomplete_dir = "/media/advseq/.run_status/pending_db_run_import"
  
version = "v1.1"
numcore = 24
padding = 0#change for ex. 150bp pad to bed regions

##turn on/off parts of this script and other options
do_clobber = True#allow overwrite of existing files
do_clean = False#remove intermediate unneeded files (!!!unimplemented!!!)
do_run = True#do run picard hsmetrics, samtools, bedtools or assume files exist
do_print = True#print run information
do_throw = True#stop execution in the case of an error.  otherwise warn

do_use_rundef = True#or run all fastqs in dir, !UNIMPLEMENTED!

do_demultiplex = False#run bcl2fastq
do_checksum = False#md5sum fastqs if bcl2fastq

do_write_qc = True#do write the final qc file if qc enabled
do_mv_output = True#move from the directory final output dir of dragen to /media/bams
do_stage_output = True#output dragen to directory to be copied later

##variant calling
do_call = False#do variant calling
do_vqsr = False#do VQSR if variant caller enabled
do_depthofcov = False#skip dragen depth of coverage or run it (unused for this qc, but makes file)

##qc 
do_qc = False#do qc process

#qc: alignment
do_hsmetrics = True#do alignment qc if qc enabled
do_run_hsmetrics = True#do alignment but assume existing hsmetrics output (do run picard runs the exec)

#qc: coverage (#bases over thresh, mean, min coverage)
do_cov = True#do coverage qc if qc enabled
do_run_cov = True#assume existing bedtools coverage -d output or run bedtools executable
do_run_cq = True

#qc: count variant call types
do_count = True#do count indels, snps in vcfs if qc enabled

#vcf concordance (picard GenotypeConcordance)
do_concordance = True#do vcf concordance qc if qc enabled
do_run_genoconcord = True
do_slicendice = True#slice by control bed file and truth bedfile to the call vcf

#RunCompletionStatus.xml (cluster density, yield, pass filter for qc)
do_xml_yield = True#get cluster density, yield, pass filter from RunCompletionStatus.xml, if processing takes <~1.5 hours after RTAComplete.txt, waits for file to generate.
rcswaitlen = 9000#seconds, 2.5 hours wait before error for RunCompletionStatus.xml to generate
rcswaitinterval = 120#two minute sleep time while waiting


##files
dbsnp_hg37 = "/staging/data/files/dbsnp_138.b37.vcf"
dbsnp_hg38 = "/staging/data/files/dbsnp_151.b38.fix.mt.vcf"
vqsr_known = ["INDEL,12.0,/staging/data/files/Mills_and_1000G_gold_standard.indels.b37.vcf",\
					"SNP,15.0,/staging/data/files/1000G_phase1.snps.high_confidence.b37.vcf", \
					"SNP,15.0,/staging/data/files/hapmap_3.3.b37.vcf", \
					"SNP,12.0,/staging/data/files/1000G_omni2.5.b37.vcf"]

#reference files
ref_hg37 = "/staging/data/ref/hg37"#dir with hash tbl for build hg19
ref_hg38 = "/staging/data/ref/hg38"#dir with hash tbl for build hg19
dragen_hg37_fa = "/staging/data/ref/hg37/human_g1k_v37.fasta"
dragen_hg38_fa = "/staging/data/ref/hg38/hg38.nochr.mt.fasta"
ref_hg37_fa = "/media/bams/workspace/advseq3.3/human_g1k_v37.fasta"
ref_hg38_fa = "/media/bams/workspace/advseq3.3/hg38.nochr.mt.fasta"
ref_genome_bed = "/staging/data/roi/hg37.bed"#to use for whole genome roi (no roi)
ref = None#reference build set in args below (hg37, hg38)


#for pass/fail % in qc coverage
qc_cut_mincov = 20
qc_cut_meancov = 30
filter_flag = 3840#samtools filter flag to remove bad reads for coverage calc

#header for qc file (see the unused below and delete if using this one)
qc_header = ["Run_Name", "Flowcell", "Well_Number", "Version", "Cluster_Density", "Pass_Filter", "Q30", "Sample_Data_Yield", "Read_Count", "Aligned_PhiX", "Aligned_Reads", "Duplicate_Rate", "Coverage_8X", "Coverage_20X", "On-Target", "Mean_Target_Coverage", "Target_Pass_Rate", "Off-Target", "Number_of_SNV", "Number_of_Indel", "Genotype_Variant_Concordance"]
#qc_header = ["Run_Name", "Flowcell", "Well_Number", "Version", "Read_Count", "Aligned_Reads", "Duplicate_Rate", "Coverage_8X", "Coverage_20X", "On-Target", "Mean_Target_Coverage", "Target_Pass_Rate", "Off-Target", "Number_of_SNV", "Number_of_Indel", "Genotype_Variant_Concordance"]


#outcol = [ None ]* len(qc_header)#qc output array for each column


##!!!!need to get runid and flowcell different way if not providing illumina seq dir  RunParameters.xml!!!
#!!!unused if qc_header used above
run_colnames = ["Run_Name", "Flowcell", "Well_Number", "Version"]
aln_colnames = 	["Read_Count", "Aligned_Reads", "Duplicate_Rate", "Coverage_8X", "Coverage_20X", "On-Target", "Mean_Target_Coverage", "Target_Pass_Rate"]
vcf_colnames = ["Number_of_SNV", "Number_of_Indel", "Genotype_Variant_Concordance"]



#locations
tmpdir = "/staging/run"#processing dir for dragen
baseseqdir = "/media/Informatics_Incoming"#bcls
rundefdir = "/media/advseq/rundef"#rundef locations
baserundir = "/media/bams"#run directory base
stagingdir = "/staging/transfer"
basepaneldir = "/media/bams/workspace/advseq3.3/"#resources dir with the panel dir (/opt/advseq/NMD0, ex)
basemovedir = "/media/bams"
bcltofq_exec = "/media/bams/workspace/advseq3.3/bin/bcl2fastq"#full path to bcl2fastq
#bcltofq_exec = "/media/bams/workspace/group_bi/dist/illumina/bcl2fastq"#full path to bcl2fastq
bedtools_exec = "/media/bams/workspace/advseq3.3/bin/bedtools"#full path to bedtools
samtools_exec = "/media/bams/workspace/advseq3.3/bin/samtools"


#key variables
ctrls = None#control wells
rundef = None#path to rundef
samplesheet = None#path to samplesheet
roi = None#path to roi
roi_concord = None#roi for vcf concordance

truth_vcf = None#vcf concordance, defined below in panel directory
truth_bed = None#vcf concordance, defined below in slice 'n dice section of vcf concordance


#file extensions used to create filenames or identify output files to use 
roi_ext = ".roi.bed"#bed file
final_vcf_ext = ".vcf"#final vcf out of dragen to use for qc (if vqsr used change to .vqsr.vcf)
run_ctrl_vcf_ext = ".vcf"#run's control sample vcf.  change to any entention with sample number
ctrl_vcf_ext = ".ctrl.vcf"#vcf for high-confidence truth vcf
ctrl_roi_ext = ".ctrl.bed"#bed for high-confidence truth intervals

#qc output filename
qc_file = None#name of outputted qc file (default runname_flowcell.qc.txt, see defining extension below)
qc_outfile_suffix = ".txt"#output file extension
qc_outfile_prefix = "QC_"#prepend to filename

#regex to identify if a file is an R1 fastq for gathering fastq pairs
re_fqgz = re.compile("^.+_R1_.+\.fastq\.gz$")



##initialize objects
ns_obj = qdx_util.Nextseq()
rundef_obj = qdx_util.Rundef()


##arguments
parser = argparse.ArgumentParser()
parser.add_argument('-seqdir', help='illumina sequence directory from: ' + baseseqdir + ' (required)')
parser.add_argument('-qc', help='run quality control report (one or more of bcl2fastq, vc, and qc)', default=False, action="store_true")
parser.add_argument('-vc', help='run variant caller.  must be on dragen machine (one or more of bcl2fastq, vc, and qc)', default=False, action="store_true")
parser.add_argument('-bcl2fastq', help='run bcl2fastq (one or more of bcl2fastq, vc, and qc)', default=False, action="store_true")
parser.add_argument('-package_fastq', help='rename fastqs for sample wells with the sample accession found in the rundef.', default=False, action="store_true")
parser.add_argument('-warn_only', help='instead of error, warn only and proceed to next sample.', default=False, action="store_true")
parser.add_argument('-flag_complete', help="create directory to flag the analysis is complete for the database (default False)", default=False, action="store_true")
parser.add_argument('-baserundir', help='final output and run directory (/media/bams/ as ex., optional, default: '+baserundir+')', default=baserundir)
parser.add_argument('-stagedir', help='directory for dragen output to be staged for copy')
parser.add_argument('-runid', help='run id wuth runname and flowcell (ex. WES011111111_Y2KK2YGXX, required if no seqdir provided)')
parser.add_argument('-panel', help='four letter panel id (ex NBDX, WES0, required if disabling the xml/runParameters.xml parse)')
parser.add_argument('-paneldir', help='directory containing panel specific files (ex.  0_based.bed, SampleSheet.csv, typically /opt/advseq/NBDX, for ex)')
parser.add_argument('-roi', help='roi bed file of regions to limit sequence to (optional, default in panel directory)')
parser.add_argument('-concordance_roi', help='roi bed file of regions to limit sequence to in the vcf concordance (optional, default is standard roi)')
parser.add_argument('-samplesheet', help='SampleSheet.csv location (optional, default in panel directory, then original sequencing directory seqdir)')
parser.add_argument('-experiment', help='run name (typically in /media/bams) with a run experiment name and the flowcell id (optional, derived from runParameters.xml by default, ex. NBDX00111111000)')
parser.add_argument('-flowcell', help='Flowcell id (nextseq w/o leading "A", optional unless norunparam flag invoked, ex. H7LLL27JX)')
parser.add_argument('-norundef', help='run on all found fastqs/bams without a rundef to indicate samples and controls (default False, uses rundef)', default=False, action="store_true")
parser.add_argument('-norunparam', help='do not use Illumina RunParameters.xml to derive the panel and flowcell id (-flowcell, -panel required then)', default=False, action="store_true")
parser.add_argument('-dowait_rcs', help='qc: wait for Illumina RunCompletionStatus.xml if not exist created approx 1.5hrs after RTAComplete.txt.  wait length: ' + str(rcswaitlen) + ' sec.', default=False, action="store_true")
parser.add_argument('-noroi', help='instead of limiting regions by bed file, utilize the entire genome(optional, default=panel_roi)', default=False, action="store_true")
parser.add_argument('-hg38', help='use genome build hg38 instead of hg37', default=False, action="store_true")
parser.add_argument('-checkwes0', help='check for panel name "WES0" to use genome build hg38 instead of hg37', default=False, action="store_true")
parser.add_argument('-notools', help='qc: assume files exist (false) or run hsmetrics, samtools, bedtools (true)', default=False, action="store_true")
parser.add_argument('-dotool_align', help='qc: assume alignment files exist (.hsmetrics) or run picard hsmetrics', default=False, action="store_true")
parser.add_argument('-dotool_coverage', help='qc: run bedtools/samtools for coverage or assume coverage files exist (.cov.bed)', default=False, action="store_true")
parser.add_argument('-dotool_cq', help='qc: run script to tally coverage output from bedtools (needs .cov.bed files) or assume cq tallied files exist (.cq.tbl)', default=False, action="store_true")
parser.add_argument('-version', help='version of this software (optional)', default=str(version))
parser.add_argument('-logfile', help='path to log file to write, default to stderr')#default to stderr
arg = parser.parse_args()

#confirm config exists
#if not os.path.exists(config):
#	message = "ERROR: configuration file specifies not found: " + config + "\n"
#	raise Exception(message)

#log file
if arg.logfile is not None:
	with open(arg.logfile) as handle:
		logfh = handle
else:
	logfh = sys.stderr

#OPTIONAL ARGS
##set reference build
##done below when checking for wes0 if specified

#set fasta to local dragen if variant call
if arg.vc:#set to dragen reference fastas
	ref_hg37_fa = dragen_hg37_fa
	ref_hg38_fa = dragen_hg38_fa


#REQUIRED ARGS
#need an action to perform, bcl2fastq, variant call (vc) or quality report (qc) or rename and md5 of sample fastqs
if not arg.qc and not arg.vc and not arg.bcl2fastq and not arg.package_fastq:#error no action chosen
	parser.print_usage()
	message = "ERROR: one or more of -qc, -vc, -bcl2fastq, -package_fastq must be chosen.\n"
	logfh.write(message)
	if logfh is not sys.stderr:
		sys.stderr.write(message)

	exit(1)#raise Exception(message)

if arg.bcl2fastq and arg.seqdir is None:
	parser.print_usage()
	message = "ERROR: -bcl2fastq requires that -seqdir is defined.\n"

	logfh.write(message)
	if logfh is not sys.stderr:
		sys.stderr.write(message)
	
	exit(1)

do_qc = arg.qc
do_call = arg.vc
do_demultiplex = arg.bcl2fastq
do_pack_fq = arg.package_fastq

#throw an error and crash on error.  warn only and continue to next sample otherwise
do_throw = not arg.warn_only

#find run name and flowcell id
#raw illumina sequence directory with bcls and xmls
if arg.seqdir is None:
	runid = None
	seqdir = None
else:
	seqdir = os.path.join(baseseqdir, arg.seqdir)
	
	if not os.path.exists(seqdir):
		parser.print_usage()
		message = "ERROR: sequence directory (-seqdir) not found: " + str(seqdir) + ", provided: " + str(arg.seqdir) + "\n"
		logfh.write(message)
		exit(1)#raise Exception(message)

	logfh.write("#parsing RunParameters.xml\n")

	rpxml = os.path.join(seqdir, "RunParameters.xml")
	if os.path.exists(rpxml):
		runname = ns_obj.experiment_name(rpxml)
		flowcell = ns_obj.flowcell_id(rpxml)
		
		if runname is None:#runname didn't parse
			message = "ERROR: couldn't gather run name (NBDX000001, etc) from RunParameters.xml: " + rpxml + "\n"
			logfh.write(message)
			raise Exception(message)

		if flowcell  is None:#flowcell id didn't parse
			message = "ERROR: couldn't gather flowcell id from RunParameters.xml: " + rpxml + "\n"
			logfh.write(message)
			raise Exception(message)

		runid = runname + "_" + flowcell
	else:
		message = "ERROR: cannot find RunParameters.xml in seq directory: " + str(seqdir) + "\n"
		logfh.write(message)
		raise Exception(message)

if arg.norundef and do_pack_fq:
	#rename of fastq with accession requires a rundef so -norundef option invalid with -accession_fastq
	message = "ERROR: rename of fastqs with rundef accession (-package_fastq) cannot be used with option (-norundef).\n"
	logfh.write(message)
	parser.print_usage()
	exit(1)

if arg.runid is not None:
	#do/not use illumina xml to get run name and flowcell id to make final directory (requires all three provided as option)
	runid = arg.runid

	runparts = runid.split("_")
	if len(runparts) != 2:
		message("ERROR: runid must consist of run name and flowcell separated by underscore, ex. WES011111111_Y2KK2YGXX\n")
		raise Exception(message)

	runname = runparts[0]
	flowcell = runparts[1]

elif arg.experiment is not None and arg.flowcell is not None:
	#requires runname and flowcell id if not looking in original xml to extract, or extracted below
	runname = arg.experiment
	flowcell = arg.flowcell
	runid = runname + "_" + flowcell


#error if run id not formed, double check found runname and flowcell
if runid is None:
	message = "ERROR: The runid was not established either from -seqdir or -runid. run name (experiment) or flowcell unknown.\n"
	
	if logfh is not sys.stderr:
		logfh.write(message)
	
	raise Exception(message)

else:

	if not os.path.exists(baserundir):#check if base processing dir exists (/media/bams, etc)
		message = "#ERROR: base run directory not found: " + baserundir + "\n"
		logfh.write(message)
		raise Exception(message)

	runpath = os.path.join(baserundir, runid)


#double check have output dir
if not os.path.exists(runpath):#make directory if doesn't exist

	if do_clobber:#allow overwrite of existing files
		message = "#WARNING: run path directory already exists.  Existing files will be overwritten: " + runpath + "\n"
		logfh.write(message)

		os.mkdir(runpath)

	else:#error if exists
		message = "#ERROR: run path directory already exists.  do_clobber set to safe.  stopping so no files overwritten: " + runpath + "\n"
		logfh.write(message)
		raise Exception(message)


####
####start process

##rundef, xmls, SampleSheet, panel id and metadata of wellnums of samples and controls
logfh.write("####\n##Collecting files and sample information.\n####\n")



##panel id and panel directory (ex. NBDX, /opt/advseq/NBDX)
panel = runname[0:4] if arg.panel is None else arg.panel
paneldir = os.path.join(basepaneldir, panel) if arg.paneldir is None else arg.paneldir


##!!!set reference build to hg38 if -hg38 or -checkwes0 and panel == WES0
if arg.hg38 or (arg.checkwes0 and panel == "WES0"):
	ref = ref_hg38
	ref_fa = ref_hg38_fa
	dbsnp = dbsnp_hg38
else:
	ref = ref_hg37
	ref_fa = ref_hg37_fa
	dbsnp = dbsnp_hg37


#use/don't use and provide an roi bed file.  if noroi it uses a genome bed
#error if roi provided and indicated to not limit by roi
if arg.noroi and arg.roi is not None:#skip roi but also provided roi.  error.
		
	message = "ERROR: conflicting arguments.  ROI file provided yet -noroi option involked. Exiting...\n"
	logfh.write(message)
	raise Exception(message)

elif arg.noroi:
	roi = "genome.roi.to.be.defined"

elif arg.roi is not None:
	roi = arg.roi#none or provided

else:
	roi = os.path.join(paneldir, panel + roi_ext)#find panel directory

if roi is None:
	message = "#ERROR: no roi determined after looking.  panel: " + str(panel) + ", panel dir: " + str(paneldir) + "\n"
	logfh.write(message)
	raise Exception(message)

#check roi exists
if not os.path.exists(roi):
	message = "#ERROR: roi file not found: " + roi + "\n"
	logfh.write(message)
	raise Exception(message)


if arg.concordance_roi is None:#vcf concordance specific roi
	roi_concord = roi
else:
	roi_concord = arg.concordance_roi

if not os.path.exists(roi_concord):
	message = "ERROR: provided concordance roi not found: " + str(roi_concord) + "\n"
	raise Exception(message)

#sample sheet, provide as arg or looks for SampleSheet.csv in the panel directory
samplesheet = os.path.join(paneldir, "SampleSheet.csv") if arg.samplesheet is None else arg.samplesheet

if do_demultiplex and not os.path.exists(samplesheet):
	message = "#ERROR: SampleSheet.csv required for bcl2fastq not found in panel (" + panel + "), file: " + samplesheet + "\n"
	logfh.write(message)
	raise Exception(message)


if arg.vc and do_stage_output:
	xferdir = os.path.join(stagingdir, runid)

	#dragen: make staging dir for transfer
	if not os.path.exists(xferdir):
		os.mkdir(xferdir)


#rundef
if not arg.norundef:
	rundef = os.path.join(rundefdir, runname + ".txt")
	if os.path.exists(rundef):
		smpls = rundef_obj.distill(rundef)
		ctrls = rundef_obj.controls#list of control well numbers

	else:
		message = "#ERROR: rundef file not found even though option to disable not set: " + rundef + "\n"
		logfh.write(message)
		raise Exception(message)
else:#DIE CUZ UNIMPLEMENTED
	rundef = None#filename
	smpls = None#dict with samples as keys (well numbers)
	ctrls = None#list

	if arg.qc or arg.vc:#need rundef for qc and variant calling (not bcl2fastq)
		message = "ERROR: option to run without a rundef (-norundef) is not implemented.  psyched you out.\n"
		logfh.write(message)
		raise Exception(message)

#dont run picard hsmetrics, samtools, bedtools, assume files exist
if arg.notools:#assume has the output files already (.hsmetrics, .cov.bed, .cq.tbl)
	do_run_hsmetrics = False
	do_run_cov = False
	do_run_cq = False

if arg.dotool_align is True:#run specific tools especially if arg notools, run alignment hsmetrics
	do_run_hsmetrics = True

if arg.dotool_coverage is True:#run samtools, bedtools coverage, cq_stats script for coverage
	do_run_cov = True
	do_run_cq = True

if arg.dotool_cq is True:
	do_run_cq = True

logfh.write("##INFO: inputs.\n")
logfh.write("#source directory: " + str(seqdir) + "\n")
logfh.write("#staging directory: " + str(stagingdir) + "\n")
logfh.write("#destination directory: " + str(runpath) + "\n")
logfh.write("#panel directory: " + str(paneldir) + "\n")
logfh.write("#rundef: " + str(rundef) + "\n")
logfh.write("#samplesheet: " + str(samplesheet) + "\n")
logfh.write("#roi: " + str(roi) + "\n")
logfh.write("#panel: " + str(panel) + "\n")

logfh.write("####\n##INFO: work to be done.\n")
logfh.write("#do use rundef file: " + str(not arg.norundef) + "\n")
logfh.write("#do limit by roi: " + str(not arg.noroi) + "\n")
logfh.write("#do bcl2fastq: " + str(do_demultiplex) + "\n")
logfh.write("#do variant calling: " + str(do_call) + "\n")
logfh.write("#--do vqsr: " + str(do_vqsr) + "\n")
logfh.write("#do qc: " + str(do_qc) + "\n")
logfh.write("#--do qc alignment: " + str(do_hsmetrics) + "\n")
logfh.write("#--do qc coverage: " + str(do_cov) + "\n")
logfh.write("#--do count variant types (indel, snp): " + str(do_count) + "\n")
logfh.write("#--do qc vcf concordance: " + str(do_concordance) + "\n")
logfh.write("#--do slice vcf by roi: " + str(do_slicendice) + "\n")
logfh.write("#do write qc to file: " + str(do_write_qc) + ".\n")

logfh.write("#do overwrite existing: " + str(do_clobber) + "\n")
logfh.write("#do clean up intermediate files: " + str(do_clean) + "\n")
logfh.write("#do run tools (not dry run) " + str(do_run) + "\n")


##bcl2fastq 
if do_demultiplex:#do demultiplexing
	logfh.write("BEGIN bcl2fastq conversion\n")

	demulti_obj = qdx_util.Demultiplex(executable=bcltofq_exec)

	res = demulti_obj.bcl2fastq(seqdir, runpath, samplesheet)

	print "bcl2fastq output\n" + str(res) + "\nendbcl2fastq output\n"


	if do_checksum:#md5sum all sample fastqs
		print "do checksum"
		fqs = list()

		dircontents = os.listdir(runpath)
		for dircontent in dircontents:
			if re_fqgz.search(dircontent) is not None:

				for wellnum in smpls:
					if dircontent.split("_")[0] == wellnum:#found sample well for read 1
						r1 = os.path.join(runpath, dircontent)
						r2 = os.path.join(runpath, dircontent.replace("_R1_", "_R2_"))

						if os.path.exists(os.path.join(runpath, r2)):#check for read 2
							fqs.append((wellnum, r1, r2))
						else:
							message = "ERROR: found read 1, but read 2 not found.\n\tread 1: " + r1 + "\n\tread 2: " + r2 + "\n"
							logfh.write(message)
							raise Exception(message)

						cmd = "md5sum " + r1
						cs = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						res = cs.communicate()

						if cs.returncode != 0:
							message = "ERROR: error occurred when running command.  returned: " + str(cs.returncode) + ", cmd: " + cmd + "\n"
							logfh.write(message)
							raise Exception(message)

						r1_md5 = os.path.join(runpath, r1 + ".md5")
						with open(r1_md5, 'w') as handle:
							handle.write(str(res[0]) + "\n")
							
						#check stderr!!!
						logfh.write("INFO: md5sum complete: " + r1_md5 + "\n")

						cmd = "md5sum " + r2
						cs = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						res = cs.communicate()

						if cs.returncode != 0:
							message = "ERROR: error occurred when running command.  returned: " + str(cs.returncode) + ", cmd: " + cmd + "\n"
							logfh.write(message)
							raise Exception(message)

						r2_md5 = os.path.join(runpath, r2 + ".md5")
						with open(r2_md5, 'w') as handle:
							handle.write(str(res[0]) + "\n")

						#check stderr!!!
						logfh.write("INFO: md5sum complete: " + r2_md5 + "\n")

	logfh.write("DONE demultiplex\n")
else:
	logfh.write("SKIP demultiplex\n")


if do_call:#do variant calling
	logfh.write("BEGIN variant calling\n")
	drgn_outdir = xferdir if do_stage_output else runpath

	if not os.path.exists(drgn_outdir):#make output dir
		os.mkdir(drgn_outdir)

	drgn = dragen.Call(call_bed=roi, reference_dir=ref, output_dir=drgn_outdir, do_vqsr=do_vqsr, vqsr_known=vqsr_known, dbsnp=dbsnp, do_clobber=do_clobber, nt=numcore, tmpfile_dir=tmpdir, do_cov=do_depthofcov, padding=padding)

	dircontents = os.listdir(runpath)

	fqs = list()#list of 3tuples of wellnum, r1, r2 fastqs to run with dragen
	for dircontent in dircontents:
		if re_fqgz.search(dircontent) is not None:

			for wellnum in smpls:
				if dircontent.split("_")[0] == wellnum:#found sample well for read 1
					r1 = dircontent
					r2 = dircontent.replace("_R1_", "_R2_")

					if os.path.exists(os.path.join(runpath, r2)):#check for read 2
						fqs.append((wellnum, r1, r2))
					else:
						message = "ERROR: found read 1, but read 2 not found.\n\tread 1: " + r1 + "\n\tread 2: " + r2 + "\n"
						logfh.write(message)
						raise Exception(message)

	
	#die unless have all sample fastqs
	if len(smpls.keys()) > len(fqs):
		message = "ERROR: did not find all fastqs for sample well numbers.\n\tsample well numbers: " + ", ".join(smpls.keys()) + "\n\tfastqs: " + str(fqs) + "\n"
		logfh.write(message)
		raise Exception(message)

	
	for readpair in fqs:
		wellnum = readpair[0]
		r1 = os.path.join(runpath, readpair[1])
		r2 = os.path.join(runpath, readpair[2])

		cmd = drgn.call_from_pair_sex(r1, r2, output_prefix=wellnum)
		print cmd
		if do_print:
			logfh.write("INFO command: " + cmd + "\n")

		if do_run:
			logfh.write("INFO VARIANTCALL sample start: " + readpair[0] + "\n")
			res = drgn.run(cmd)
			#print str("".join(res))
			logfh.write("INFO VARIANTCALL sample complete: " + readpair[0] + "\n")
		else:
			message = "INFO: do_run switch indicates this is a dry run and will not be executed on dragen: " + readpair[0] + "\n"
			logfh.write(message)
		
	if do_mv_output:#move all output files to their final resting place
		#destination = os.path.join(basemovedir, panel)

		message = "INFO: moving dragen output files...\nfrom: " + str(stagingdir) + "\n\tto: " + str(runpath) + "\n"
		logfh.write(message)

		if not os.path.exists(drgn_outdir):
			message = "ERROR: directory of final dragen output files to move not found: " + drgn_outdir + "\n"
			raise Exception(message)

		dircontents = os.listdir(drgn_outdir)
		if len(dircontents) == 0:
			message = "ERROR: no files in staging directory to move off dragen to final run directory.  from: " + str(drgn_outdir) + ", to: " + str(runpath) + "\n"
			logfh.write(message)
			raise Exception(message)

		#for dircontent in dircontents:
		try:
			logfh.write("Copying dragen output from source to destination: src: " + drgn_outdir + ", dest: " + runpath + "\n")
				#fpath = os.path.join(drgn_outdir, dircontent)
				#print "to copy but is not: " + os.path.join(fpath, dircontent) + ", to: " + destination
				#shutil.copy(os.path.join(fpath, dircontent), destination)
				#logfh.write("Success.  Moved dragen output from source to destination: src: " + fpath + ", dest: " + runpath + "\n")

			#add switch '-p' to maintain timestamp on filecopy
			cmd = "cp -rf " + drgn_outdir + "/* " + runpath
			sp = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
			sp.wait()

			if sp.returncode != 0:
				message = "ERROR: copying from staging dir to final dir had non-zero error code: " + str(sp.returncode) + ", staging: " + drgn_outdir + ", final: " + runpath + "\n"
				raise Exception(message)

			cmd = "rm -r " + drgn_outdir
			sp = subprocess.Popen(cmd, stdout=sys.stderr, stderr=sys.stderr, shell=True, preexec_fn=os.setpgrp)
				
		except:
			message = "ERROR: cannot move dragen output file: " + drgn_outdir + "\n"
			logfh.write(message)
			raise Exception(message)
		
	else:
		message = "INFO: No move from staging directory indicated from: " + str(drgn_outdir) + "\n"
		logfh.write(message)

if do_call:
	message = "INFO: final location of output files: " + str(runpath) + "\n"
	logfh.write("DONE variant calling\n")

else:
	logfh.write("SKIP variant calling\n")


if do_qc:##run qc report
	logfh.write("BEGIN QC\n")

	smpl_aln = list()
	#smplout = dict()#outputted columns
	smpl_order = list()#keep order of samples in the bams
	colnames = list()
	gconcord_avg = None

	if ctrls is None:
		message = "ERROR: no controls found upon starting qc.  rundef may not have a control or it did not parse correctly.  rundef: " + rundef + "\n"
		logfh.write(message)
		
		if logfh is not sys.stderr:
			sys.stderr.write(message)
		#raise or no?#raise Exception(message)

	####
	##find each bam for all samples (including controls)
	bams = list()
	for smpl, smplinfo in smpls.iteritems():
		#add sample to have an entry in output
		#smplout[smpl] = list()#list of columns
		smpl_order.append(smpl)

		#look for bams, make optional for vcf concordance!!!
		if os.path.exists(os.path.join(runpath, smpl+".bam")):
			bams.append(os.path.join(runpath, smpl+".bam"))
		else:
			message = "ERROR: expected bam file not found: " + os.path.join(runpath, smpl+".bam") + "\n"
			logfh.write(message)
			raise Exception(message)

	if len(bams) < len(smpls.keys()) and (do_hsmetrics or do_cov):
		message = "ERROR: not all bams found in sample rundef gathered.\n\tsamples: " + str(smpls.keys()) + "\n\t" + str(bams) + "\n"
		logfh.write(message)
		raise Exception(message)

	####
	##alignment qc
	if do_hsmetrics:
		logfh.write("BEGIN ALIGNMENT QC\n")
		picard_aln_obj = picard_util.AlignUtil(do_run=do_run_hsmetrics)

		for bam in bams:
			#no reference used since no dropout evaluated#colvals = picard_aln_obj.qc_report_stats(bam, roi, reference=ref_fa, do_clean=do_clean)
			colvals = picard_aln_obj.qc_report_stats(bam, roi, do_clean=do_clean)
			smpl_aln.append(colvals)

			aln_colnames = picard_aln_obj.colnames

		logfh.write("DONE ALIGNMENT QC\n")

	####
	##coverage qc, pass filter
	if do_cov:
		logfh.write("BEGIN COVERAGE QC\n")

		cqs = list()#list of values with %failing coverage 

		for bam in bams:#change from sp of script to pull in of method!!!
			bt_cov_file = bam + ".cov.bed"
			cq_outfile = bam + ".cq.tbl"

			#run bedtools coverage or assume output files already exist (skip)
			if do_run_cov:
			
				#filter bad alignments with filter flag (3844, or above)
				cmd_st = samtools_exec + " view -u -F " + str(filter_flag) + " " + bam 

				#run bedtools coverage -d
				cmd_bt = bedtools_exec + " coverage -d -abam stdin -b " + roi
	

				if do_print:
					print cmd_st + " | " + cmd_bt + " > " + bt_cov_file

				with open(bt_cov_file, 'w') as handle:
					#samtools rm dup, 0map, etc reads for input into bedtools
					st = subprocess.Popen(cmd_st, stdout=subprocess.PIPE, stderr=sys.stderr, bufsize=1, shell=True)
					bt = subprocess.Popen(cmd_bt, stdin=st.stdout, stdout=handle, stderr=sys.stderr, shell=True)
					st.stdout.close()
					bt.communicate()

				#if res is not None and res[1] != "":#res[1]
				if bt.returncode != 0: 
					message = "ERROR: response from samtools with bedtools coverage was an exception. status: " + str(bt.returncode) + "\n"
					#for line in res[1]:
					#	message += "ERROR: " + line + "\n"
					#message += "ERROR: " + str(res) + "\n"
					logfh.write(message)
					raise Exception(message)

			if do_run_cq:#make cq coverage stats or existing already
				#make cq stats table (min, meancov, total, >20X)
				cmd = "cq_cov_qc.py " + bt_cov_file + " " + roi

				if do_print:
					print cmd

				with open(cq_outfile, 'w') as handle:
					cq = subprocess.Popen(cmd, stdout=handle, stderr=subprocess.PIPE, shell=True, bufsize=1)
					res = cq.communicate()

				#if res[1] is not None and res[1] != "":
				if cq.returncode != 0:
					message = "ERROR: response from cq_cov_stat.py was an exception. exit status: " + str(cq.returncode) + "\n"
					message += "ERROR: " + str(res) + "\n"
					logfh.write(message)
					raise Exception(message)

			#with open(cq_outfile) as handle:
			#parse cq stats
			with open(cq_outfile) as handle:
				lines = handle.readlines()
			
			failcov = 0
			numregion = 0
			over8 = 0
			over20 = 0
			totalbases = 0
			meancov_total = 0
			mincov_total = 0

			#currname + "\t" + str(cov_avg_form) + "\t" + str(cov_min) + "\t" + str(cov_under2) + "\t" + str(cov_under) + "\t" + str(cov_len)
			for line in lines:
				cols = line.rstrip().split("\t")

				meancov = float(cols[1])#per region
				meancov_total += meancov#sum of mean cov over each region

				mincov = float(cols[2])
				mincov_total += mincov

				#under8x = int(cols[3])/float(totalbases)
				#under20x = int(cols[4])/float(totalbases)
				#if int(cols[3]) > 8:
				over8 += int(cols[3])#num bases > 8X
				#if int(cols[4]) > 20:
				over20 += int(cols[4])#num bases > 20X

				totalbases += int(cols[5])#length of region

				if mincov < qc_cut_mincov or meancov < qc_cut_meancov:
					failcov += 1

				numregion += 1

			if numregion > 0:
				failcov_avg = float(failcov)/numregion
				passcov_avg = 1-failcov_avg

				meancov_avg = float(meancov_total)/numregion
				mincov_avg = float(mincov_total)/numregion
			else:
				message = "ERROR: number of regions found in qc coverage is 0. file: " + str(cq_outfile) + "\n"
				logfh.write(message)
				raise Exception(message)

			#for qc report %over 8x,20x, and pass fail % based on mean min cov
			#outcol[7] = str(over8/float(numregion))
			#outcol[8] = str(over20/float(numregion))
			#outcol[11] = str(passcov_avg)

			meancov_avgperc = round(float(meancov_avg)*100, 3)
			mincov_avgperc = round(float(mincov_avg)*100, 3)

			passcov_avgperc = round(float(passcov_avg)*100, 3)
			#over8perc = round(over8/float(numregion)*100, 3)#true percent rounded to 3 decimals
			#over20perc = round(over20/float(numregion)*100, 3)
			over8perc = round(over8/float(totalbases)*100, 3)#avg #bases over 8x divided by total bases in roi
			over20perc = round(over20/float(totalbases)*100, 3)

			form_cov = str(passcov_avgperc), str(over8perc), str(over20perc), str(meancov), str(mincov)
			cqs.append(form_cov)#turn into % passed


		logfh.write("DONE COVERAGE QC\n")


	####
	##vcf count snp, indels
	if do_count:
		logfh.write("BEGIN COUNT VARIANT TYPE\n")

		counts = list()#list of two tuples, (num snp, num indel) in smpl_order order

		for smpl in smpl_order:
			vcf = smpl + final_vcf_ext
			vcffile = os.path.join(runpath, vcf)

			if not os.path.exists(vcffile):
				if os.path.exists(vcffile+".gz"):#gunzip vcf.gz if needed
					logfh.write("INFO: found gzipped vcf, unzipping...   " + vcffile + ".gz\n")
					res = subprocess.Popen("gunzip "+vcffile+".gz", shell=True).communicate()

					if res[1] is not None:
						message = "ERROR: couldn't gunzip run control vcf: " + vcffile + "\nERROR: " + str(res[1]) + "\n"
						raise Exception(message)
				else:
					message = "ERROR: expected vcf file for sample not found. sample: " + str(smpl) + ", file: " + vcffile + "\n"
					logfh.write(message)
					raise Exception(message)

			numsnp = 0
			numindel = 0

			with open(vcffile) as handle:
				lines = handle.readlines()

			for line in lines:
				if not line.startswith("#"):
					cols = line.rstrip().split("\t")

					ref = cols[3]
					alt = cols[4]

					if alt == ".":
						continue

					if len(ref) == len(alt):
						if len(ref) == 1:#snp
							numsnp += 1

						else:#delins, currently snp
							numsnp += 1

					else:
						numindel += 1

			counts.append((numsnp, numindel))

		logfh.write("DONE COUNT VARIANT TYPE\n")

	####
	##vcf concordance qc
	if do_concordance:
		logfh.write("BEGIN VCF CONCORDANCE\n")

		if len(ctrls) == 0:#no control wells
			gconcord_avg = str()

		else:#do concordance
			##vcf file of run control well
			file_ctrl_vcf = ctrls[0] + run_ctrl_vcf_ext
			run_ctrl_vcf = os.path.join(runpath, file_ctrl_vcf)

			if not os.path.exists(run_ctrl_vcf):

				if os.path.exists(run_ctrl_vcf+".gz"):#gunzip vcf.gz if needed
					logfh.write("INFO: found gzipped vcf, unzipping...   " + run_ctrl_vcf + ".gz\n")
					sp = subprocess.Popen("gunzip "+run_ctrl_vcf+".gz", shell=True)
					res = sp.communicate()
			
					if sp.returncode != 0:#error
						message = "ERROR: could not gunzip run control vcf: " + run_ctrl_vcf
						if res[1] is None:
							message += "\n"
						else:
							message +=  ", error string: " + ref[1] + "\n"
							raise Exception(message)


				else:#no vcf found
					message = "ERROR: vcf file of the control sample not found.  sample: " + str(ctrls[0]) + ", vcf: " + run_ctrl_vcf + "\n"
					logfh.write(message)
					raise Exception(message)


			##high confidence truth vcf 
			if truth_vcf is None:#defined above or no
				truth_vcf = os.path.join(paneldir, panel + ctrl_vcf_ext)

			if not os.path.exists(truth_vcf):
				message = "ERROR: control vcf for concordance not found: " + truth_vcf + "\n"
				logfh.write(message)
				raise Exception(message)


			##intersect both truth bed file with call vcf and panel roi bed file with truth vcf
			if do_slicendice:

				#slice our calls with high confidence truth bed file
				#get truth vcf's bed interval file in panel dir
				if truth_bed is None:
					truth_bed = os.path.join(paneldir, panel + ctrl_roi_ext)

				if not os.path.exists(truth_bed):
					message = "ERROR: cannot find bed interval file for truth vcf in vcf concordance slice 'n dice: " + str(truth_bed) + "\n"
					logfh.write(message)
					raise Exception(message)

				thecall = run_ctrl_vcf + ".slice.vcf"
			
				slice1cmd = bedtools_exec + " intersect -header -wa -a " + run_ctrl_vcf + " -b " + truth_bed

				if do_print:
					print slice1cmd

				with open(thecall, 'w') as handle:
					sp = subprocess.Popen(slice1cmd, stdout=handle, stderr=subprocess.PIPE, shell=True)
					res = sp.communicate()
			
					if sp.returncode != 0:#error
						message = "ERROR: could not slice the call vcf with high-confidence roi bedfile."
						if res[1] is None:
							message += "\n"
						else:
							message +=  ", error string: " + ref[1] + "\n"
						raise Exception(message)

				#if res[1] is not None:
				#	message = "ERROR: could not slice the our call vcf with high confidence truth vcf bed.\n" + res[1] + "\n"
				#	raise Exception(message)

				#slice high confidence truth vcf with our panel roi bed
				thetruth = os.path.join(runpath, os.path.basename(truth_vcf)) + ".slice.vcf"

				slice2cmd = bedtools_exec + " intersect -header -wa -a " + truth_vcf + " -b " + roi_concord

				if do_print:
					print slice2cmd

				with open(thetruth, 'w') as handle:
					sp = subprocess.Popen(slice2cmd, stdout=handle, stderr=subprocess.PIPE, shell=True)
					res = sp.communicate()
			
					if sp.returncode != 0:#error
						message = "ERROR: could not slice the truth vcf with our roi bedfile."
						if res[1] is None:
							message += "\n"
						else:
							message +=  ", error string: " + ref[1] + "\n"
						raise Exception(message)

			else:
				thecall = run_ctrl_vcf
				thetruth = truth_vcf

			#check for call and truth vcf, sliced or not 
			if not os.path.exists(thecall):
				message = "ERROR: the call vcf (sliced) to be used for concordance not found: " + thecall + "\n"
				logfh.write(message)
				raise Exception(message)
			if not os.path.exists(thetruth):
				message = "ERROR: the truth (sliced high confidence) vcf to be used for concordance not found: " + thetruth + "\n"
				logfh.write(message)
				raise Exception(message)

			#lib
			picard_vcf_obj = picard_util.VcfConcordance(do_clean=do_clean, do_run=do_run_genoconcord)

			#detail, summary and concordance stats
			s, d, c = picard_vcf_obj.genoconcord(thecall, thetruth, outprefix=file_ctrl_vcf, do_output_vcf=False)

			#collect vals for each column separating them by row.
			#d, s, c all are arrs of two-tuples.  snp in elem 0, indel 1 for each col
			#summary_metrics snp and indel array
			s_snp = list()
			s_indel = list()

			if runname is None or flowcell is None:
				message = "ERROR: vcf concordance finds either runname or flowcell id not found. runname: " + str(runname) + ", flowcell: " + str(flowcell) + "\n"
				logfh.write(message)
				raise Exception(message)

			wellnum = picard_vcf_obj.call#sample id


			##genotype concordance average
			gconcord_snp, gconcord_indel = picard_vcf_obj.get_summary_gconcord(s)
			#gconcord_avg = (float(gconcord_snp) + float(gconcord_indel))/2
			numsnp = picard_vcf_obj.numsnp
			numindel = picard_vcf_obj.numindel
			numtotal = numsnp + numindel

	
			#average the contribution of snps to snp concordance and contribution of indel concordance
			# = gconcord of snp * total number of snps plus gconcord of indel * total indels. divide by total variants (snp+indel)
			gconcord_avg = float((float(gconcord_snp)*numsnp + float(gconcord_indel)*numindel)/numtotal)
			gconcord_avgperc = round(float(gconcord_avg)*100, 3)
			vcfcol = [str(numsnp), str(numindel), str(gconcord_avgperc)]

		logfh.write("DONE VCF CONCORDANCE\n")


	if do_xml_yield:
		#get yields cluster density, data yield, passing filter from illumina xml
		logfh.write("BEGIN XML SAMPLE YIELD\n")

		if seqdir is None:
			message = "qc specifies to get yields from illumina RunCompletionStatus.xml but no seqdir provided (-seqdir)\n"
			
			if do_throw:
				raise Exception("ERROR: " + message)
			else:
				logfh.write("WARN: " +  message + "WARN: continuing without Illumina yields\n")
				if logfh is not sys.stderr:
					sys.stderr.write("WARN: " + message)
		
			cdensity = str()
			pfilter = str()
			estyield = str()


		else:#get stats
			rcsxml = os.path.join(seqdir, "RunCompletionStatus.xml")

	                #fix for if RTAComplete.txt and waiting for RunCompletionStatus.xml to generate 1.5hr later.
	                #waits up to rcswaitlen (>1.5 hours) until fails if arg.dowait_rcs==True
	                if not os.path.exists(rcsxml) and arg.dowait_rcs is True:
	                        waittime = datetime.datetime.now() + datetime.timedelta(seconds=rcswaitlen)
	
				while datetime.datetime.now() < waittime:
					message = "WARN: waiting for RunCompletionStatus.xml to generate from sequencer (~1.5hrs after sequencing completes RTAComplete.txt).\n\tsleeping...\n"
					logfh.write(message)
					time.sleep(rcswaitinterval)
	
					if os.path.exists(rcsxml):#stop waiting if found
						break

			if not os.path.exists(rcsxml):
				message = "ERROR: RunCompletionStatus.xml not found: " + rcsxml + "\n"
				raise Exception(message)

			cdensity, pfilter, estyield = ns_obj.final_yields(rcsxml)

			if cdensity is None:
				cdensity = str()
			else:
				cdensity = str(round(float(cdensity), 3))

			if pfilter is None:
				pfilter = str()
			else:
				pfilter = str(round(float(pfilter), 3))

			if estyield is None:
				estyield = str()
			else:
				estyield = str(round(float(estyield), 3))
			
		logfh.write("DONE XML SAMPLE YIELD\n")

	else:
		cdensity = str()
		pfilter = str()
		estyield = str()

	##qc output columns from header
	#1#Run_Name
	#2#Flowcell
	#3#Well_Number
	#4#Version
	#5#Cluster_Density
	#6#Pass_Filter
	#7#Q30
	#8#Sample_Data_Yield
	#9#Read_Count
	#10#Aligned_PhiX
	#11#Aligned_Reads
	#12#Duplicate_Rate
	#13#Coverage_8X
	#14#Coverage_20X
	#15#On-Target
	#16#Mean_Target_Coverage
	#17#Target_Pass_Rate
	#18#Off-Target
	#19#Number_of_SNV
	#20#Number_of_Indel
	#21#Genotype_Variant_Concordance

	#qc output, merge into one
	outlines = list()
	for i, smpl in enumerate(smpl_order):
		outcol = [ None ]* len(qc_header)#qc output array for each column

		#outcol = [runname, flowcell, smpl, version]#insert padded cols for Cluster_Density, Pass_Filter", "Q30", "Sample_Data_Yield"
		outcol[0] = runname
		outcol[1] = flowcell
		outcol[2] = smpl
		outcol[3] = version
		outcol[4] = cdensity#cluster_density
		outcol[5] = pfilter#pass_filter
		outcol[6] = str()#Q30#unfilled
		outcol[7] = estyield#Sample_Data_Yield
		outcol[9] = str()#Aligned_PhiX#unfilled

		if do_hsmetrics:

			#hs metrics array order
			#total reads=col[0], aligned reads=col[1], dup reads=col[2],on target=col3, mean cov=col4, off target=col5

			outcol[8] = smpl_aln[i][0]#read count
			outcol[10] = smpl_aln[i][1]#aligned reads
			outcol[11] = smpl_aln[i][2]#duplicate reads
			outcol[14] = smpl_aln[i][3]#on-target
			outcol[15] = smpl_aln[i][4]#mean coverage
			outcol[17] = smpl_aln[i][5]#off-target


		if do_cov:
			#cq output cols
			#str(passcov_avgperc), str(over8perc), str(over20perc), str(meancov), str(mincov)
			outcol[12] = str(cqs[i][1])#>8x
			outcol[13] = str(cqs[i][2])#>20x
			#outcol[15] = str(cqs[i][3])#mean cov
			outcol[16] = str(cqs[i][0])#pass rate

		if do_count:
			#outcol += counts[i]
			outcol[18] = counts[i][0]
			outcol[19] = counts[i][1]

		if do_concordance:

			outcol[20] = str(round(float(gconcord_avg)*100,3)) if gconcord_avg != "" else gconcord_avg


		strcol = [ str(val) for val in outcol ]#coerce all to string

		outlines.append("\t".join(strcol))


		
	#make qc header
	'''colnames = run_colnames

	#insert padded cols for Cluster_Density, Pass_Filter", "Q30", "Sample_Data_Yield", "Aligned_PhiX"
	#colnames += ["Cluster_Density", "Pass_Filter", "Q30", "Sample_Data_Yield"]

	if do_hsmetrics:
		#aln_colnames.insert(1,"Aligned_PhiX")#insert for PhiX
		colnames += aln_colnames

	if do_concordance:#add column names
		colnames += vcf_colnames

		#insert padded cols for Cluster_Density, Pass_Filter", "Q30", "Sample_Data_Yield", "Aligned_PhiX"
		#outcol += ["Cluster_Density", "Pass_Filter", "Q30", "Sample_Data_Yield"]'''

	if do_write_qc:

		if qc_file is None:
			qc_file = str()
			
			if runpath is not None:
				qc_file = runpath + "/"
			
			qc_file += qc_outfile_prefix + runname + "_" + flowcell + qc_outfile_suffix

		with open(qc_file, 'w') as handle:
			handle.write("\t".join(qc_header) + "\n")

			for line in outlines:
				handle.write(line + "\n")

		logfh.write("WROTE QC FILE: " + qc_file + "\n")

	logfh.write("DONE QC\n")
else:
	logfh.write("SKIP QC\n")


if do_pack_fq:#rename fastqs to include the accession found in rundef, md5sum the files
	logfh.write("BEGIN FASTQ RENAME and MD5 CHECKSUM\n")

	if rundef is None:#check rundef defined
		message = "ERROR: no rundef found for rename of fastqs\n"
		logfh.write(message)
		raise Exception(message)

	#sys call to script for rename and md5
	cmd = "vrntx_export.py " + rundef + " " + runpath

	renamer = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	res = renamer.communicate()

	if renamer.returncode != 0:
		message = "ERROR: error occurred when running command.  returned: " + str(renamer.returncode) + ", cmd: " + cmd + "\n"
		
		if res[1] is not None and res[1] != "":
			message += res[1].rstrip() + "\n"

		logfh.write(message)
		raise Exception(message)

	logfh.write("DONE RENAME AND CHECKSUM\n")
else:
	logfh.write("SKIP FASTQ ACCESSION RENAME\n")


if arg.flag_complete:
	logfh.write("BEGIN FLAG ANALYSIS COMPLETE.\nINFO: making directory flag: " + iscomplete_dir + "\n")
	
	donedir = os.path.join(iscomplete_dir, runid)

	try:
		os.mkdir(donedir)

		logfh.write("INFO: successfully created flagging complete. mkdir " + donedir + "\n")

	except OSError as err:
		message = "ERROR: could not create directory to flag complete: " + donedir + "\n"
		message += "ERROR: exit status " + str(err.errno) + "\n"
		message += "ERROR: error string: " + str(err.strerror) + "\n"
		logfh.write(message)

		if logfh is not sys.stderr:
			sys.stderr.write(message)

	logfh.write("DONE FLAG COMPLETE DIRECTORY CREATION\n")
else:
	logfh.write("SKIP FLAG ANALYSIS COMPLETE\n")


#parse config
#after all edico done, run hsmetrics, coverage, vcfconcordance
#merge all into qc report


#IF CONFIG IS WES0
#md5sum the fastqs
#rename bam
#set up copy



logfh.write("DONE.\n")

exit()
