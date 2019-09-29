package main

import "fmt"
import "strings"
import "flag"
import "log"
import "os"
import "time"
import "runtime"
import "io"
import "io/ioutil"
import "path/filepath"
import "bufio"
import "math/rand"
import "compress/gzip"

import "genomics/common"
import "genomics/storage"

type AlignMetadata struct {
	real_name string
	rg_id     string
	rg_pu     string
}

// Values are filled in at compilation by flags in the makefile.
var Version = "__placeholder__"
var Build = "__placeholder__"

var ngm_bin string
var trimmomatic_bin string
var picard_bin string
var bamutils_bin string
var sambamba_bin string

var nthreads int
var ref_genome string

var sambam bool

var trim_prefix string
var cull_bed string

var sample string
var sample_metadata AlignMetadata
var use_metadata bool

var ramdisk string
var use_ramdisk bool

type DuplicateOp func(string, string, int, string) ([]string, error)

func showVersion(is_log bool) {
	if is_log {
		log.Printf("align %s build %s", Version, Build)
	} else {
		fmt.Println("align " + Version + " build " + Build)
	}
}

func usage() {
	fmt.Println(`Usage: align -r REF -i JOB [opts]*

  Opts:
	 -n 	 Thread Count
	 -r REF	 Reference Genome
	 -i JOB	 Input Job File
	 -di 	 Delete Input Files
	 -lossy	 Use lossy pipeline (default false)
	 -sb	 Use sambamba tools (default false)
	 -m	 Extract metadata from file (default false)
	 -u USR	 Username
	 -p PWD	 Password
	 -creds CRDS	 Credentials file (yaml)
	 		 (username:, password:)
	 -P 	 Project Namespace (GCP)
	 -x 	 Scan Prefix (for Trimming)
	 -c FIL	 Cull Non-repetive Regions
	 -d 1/0	 Mark PCR Duplicates (default true)
	 -z 	 Output Compression
	 -w 	 Working directory (used for temporary files and downloads)
	 -o 	 Output Path
	 -ramdisk 	 Use Ramdisk (default false)
	 -l 	 Log to File (default stdout)
	 -V 	 show program version
	 -----------------------------------
	 -ngm PATH	 Path to ngm binary
	 -trimmomatic PATH	 Path to trimmomatic binary
	 -picard PATH	 Path to picard binary
	 -bamutils PATH	 Path to bamutils binary
	 -sambamba PATH	 Path to sambamba binary`)
}

// Returns (sampleName, fileurl1, hashurl1, fileurl2, hashurl2)
func parseJob(path string, creds common.AuthenticationCreds) (string, string, string, string, string, error) {

	var data map[string]common.AlignSample
	var v common.AlignSample
	var file string
	var err error

	if file, err = storage.Get(path, creds, false); err != nil {
		return "", "", "", "", "", err
	}

	if data, err = common.LoadAlignJob(file); err != nil {
		return "", "", "", "", "", err
	}

	if len(data) != 1 {
		return "", "", "", "", "", fmt.Errorf("too many jobs: %d", len(data))
	}

	for _, v = range data {
		break
	}

	return v.Name, v.Files[0][0], v.Files[0][1], v.Files[1][0], v.Files[1][1], nil
}

func initializeMetadata(mdata *AlignMetadata, s string) {
	mdata.real_name = s
	mdata.rg_id = s
	mdata.rg_pu = s
}

func extractMetadata(s string) (AlignMetadata, error) {

	var reader io.Reader
	var mdata AlignMetadata

	initializeMetadata(&mdata, sample)
	mdata.real_name = common.GetLastSubstring(sample, ".")

	fd, err := os.Open(s)
	if err != nil {
		return mdata, err
	}
	defer fd.Close()

	// Read first line from the file
	reader = fd
	if strings.HasSuffix(s, ".gz") {
		gzreader, err := gzip.NewReader(fd)
		if err != nil {
			return mdata, err
		}
		defer gzreader.Close()
		reader = gzreader
	}

	buffered := bufio.NewReader(reader)
	first, err := buffered.ReadString('\n')
	if err != nil {
		return mdata, err
	}
	first = strings.TrimRight(first, "\n")

	fparts := strings.Split(first, ":")
	if len(fparts) >= 4 {
		flowcell_id := fparts[2]
		flowcell_lane := fparts[3]
		index_seq := fparts[len(fparts)-1]

		mdata.rg_id = fmt.Sprintf("%s.%s.%s", flowcell_id, flowcell_lane, index_seq)
		mdata.rg_pu = mdata.rg_id
	}

	return mdata, nil
}

// syncs, and then saves or deletes the temporary log file
func closeLogFile(logFile string, logFd *os.File, saveLogs *string, creds *common.AuthenticationCreds) {
	if logFile != "" {
		logFd.Sync()
		logFd.Close()

		if *saveLogs != "" {
			storage.Put(*saveLogs, logFile, *creds, false)
		} else {
			os.Remove(logFile)
		}
	}
}

func cleanupTemp(tmpDir string) {
	if tmpDir != "" {
		log.Printf("Cleaning temporary files...")
		err := os.RemoveAll(tmpDir)
		if err != nil {
			log.Printf("Error cleaning up: %s", err.Error())
		}
	}
}

type multiString []string

func (i *multiString) String() string {
	return "list of strings"
}

func (i *multiString) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func dump_flag(f *flag.Flag) {
	if f.Name == "cas" {
		log.Printf("    %-30s: %s = %v (default: [])", f.Usage, f.Name, storage.CasCaches)
	} else {
		val := f.Value.String()
		defVal := f.DefValue
		if val == "" {
			val = "\"\""
		}
		if defVal == "" {
			defVal = "\"\""
		}
		log.Printf("    %-30s: %s = %s (default: %s)", f.Usage, f.Name, val, defVal)
	}
}

func mainWithCode() int {

	var ref_path string

	var creds common.AuthenticationCreds
	var user, pass, credsFile, project string

	var jobfile string
	var lossy bool

	var compression int
	var op_path, working_dir string

	var mark_repetitive bool
	var markrep_param int

	var to_log bool
	var show_version bool
	var log_file string
	var logFd *os.File

	var r1, r1h, r2, r2h string
	var s1, s2 string
	var remove_files bool
	var d1, d2 string

	var runPipe PipelineOp

	var casPaths = multiString{}
	var contentCaches = []common.Cas{}

	saveLogs := ""
	var res_files []string

	var ts_start, ref_start time.Time

	var err error

	flag.IntVar(&nthreads, "n", runtime.NumCPU(), "Thread Count")
	flag.StringVar(&ref_path, "r", "", "Reference Genome")

	flag.StringVar(&jobfile, "i", "", "Input Job File")
	flag.BoolVar(&remove_files, "di", false, "Delete Input Files")
	flag.BoolVar(&lossy, "lossy", false, "Use lossy pipeline")
	flag.BoolVar(&sambam, "sb", false, "Use sambamba tools")

	flag.BoolVar(&use_metadata, "m", false, "Use metadata")

	flag.StringVar(&user, "u", "", "Username")
	flag.StringVar(&pass, "p", "", "Password/Key")
	flag.StringVar(&credsFile, "creds", "", "Path to credentials file (yaml)")
	flag.StringVar(&project, "P", "", "Project Namespace (GCP)")

	flag.StringVar(&trim_prefix, "x", "", "Scan Prefix for Trimming")
	flag.StringVar(&cull_bed, "c", "", "Cull Non-repetitive Regions")

	flag.IntVar(&markrep_param, "d", 1, "Mark PCR Duplicates")

	flag.IntVar(&compression, "z", 5, "Output Compression")
	flag.StringVar(&op_path, "o", "", "Output Path")
	flag.StringVar(&working_dir, "w", "", "Working directory")
	flag.BoolVar(&use_ramdisk, "ramdisk", false, "Use Ramdisk")
	flag.BoolVar(&to_log, "l", false, "Log to File")
	flag.BoolVar(&show_version, "V", false, "Show program version")
	flag.StringVar(&ngm_bin, "ngm", "/home/ngm/ngm", "Path to ngm binary")
	flag.StringVar(&trimmomatic_bin, "trimmomatic", "/home/trimmomatic/trimmomatic-0.36.jar", "Path to trimmomatic binary")
	flag.StringVar(&picard_bin, "picard", "/home/picard-2.9.3.jar", "Path to picard binary")
	flag.StringVar(&bamutils_bin, "bamutils", "/home/bam", "Path to bamUtils binary")
	flag.StringVar(&sambamba_bin, "sambamba", "/home/sambamba_v0.6.6", "Path to bamUtils binary")
	flag.Var(&casPaths, "cas", "use caspath as a readonly content cache")
	flag.Parse()

	rand.Seed(time.Now().UnixNano())

	if to_log == true {
		log_file = filepath.Join(working_dir, ".align-log-"+fmt.Sprintf("%s", time.Now().Unix()))
		if logFd, err = os.Create(log_file); err != nil {
			panic(err)
		}
		log.SetOutput(logFd)
	} else {
		log_file = ""
		logFd = nil
	}
	defer closeLogFile(log_file, logFd, &saveLogs, &creds)

	if show_version {
		showVersion(false)
		return 0
	}

	mark_repetitive = true
	if markrep_param == 0 {
		mark_repetitive = false
	}

	showVersion(true)
	log.Printf("Args: %v", os.Args[1:])
	if ref_path == "" {
		log.Printf("Missing reference genome\n")
		usage()
		return 1
	}

	if jobfile == "" {
		log.Printf("Missing job file\n")
		usage()
		return 1
	}

	creds, err = common.ParseCreds(user, pass, project, credsFile)
	if err != nil {
		return 1
	}

	// Setup working directory to store all files. This will store
	// downloads as well.
	var info os.FileInfo
	info, err = os.Stat(working_dir)
	if err != nil || !info.IsDir() {
		log.Printf("invalid working directory: %s", working_dir)
		return 1
	}
	tmpDir, err := ioutil.TempDir(working_dir, ".align-")
	if err != nil {
		log.Printf("could not create temp directory inside %s: %s", working_dir, err.Error())
		return 1
	}
	storage.LocalDir = tmpDir
	log.Printf("Downloads and temporary files placed in: %s", tmpDir)

	defer cleanupTemp(tmpDir)

	// Load cache directories
	for _, casPath := range casPaths {
		cache := common.DirCas(casPath)
		if cache == nil {
			fmt.Fprintf(os.Stderr, "Could not use cache dir: %s\n", casPath)
			return 1
		} else {
			log.Printf("Using cache directory %s", casPath)
		}
		contentCaches = append(contentCaches, cache)
	}
	storage.CasCaches = contentCaches


	// ARGUMENTS PARSED

	log.Printf("Running with config:")
	flag.VisitAll(dump_flag)

	// Check job file
	sample, r1, r1h, r2, r2h, err = parseJob(jobfile, creds)
	if err != nil {
		log.Printf("Can't parse job file: %s\n", err.Error())
		usage()
		return 1
	}

	if r1 == "" || r2 == "" {
		log.Printf("Missing sample files\n")
		usage()
		return 1
	}

	// Acquire fastq files
	ts_start = time.Now()
	s1, err = getSample(r1, r1h, creds, false)
	if err != nil {
		log.Printf("Error getting the R1 sample: %s\n", err.Error())
		return 1
	}
	s2, err = getSample(r2, r2h, creds, false)
	if err != nil {
		log.Printf("Error getting the R2 sample: %s\n", err.Error())
		return 1
	}
	log.Printf("Time to get samples: %s\n", time.Since(ts_start))

	// Acquire reference files
	ref_start = time.Now()
	if ref_genome, err = getReference(ref_path, creds); err != nil {
		log.Printf("Error getting the reference: %s\n", err.Error())
		return 1
	}
	if trim_prefix != "" {
		if trim_prefix, err = storage.Get(trim_prefix, creds, false); err != nil {
			log.Printf("Error getting the trimmomatic prefix: %s\n", err.Error())
			return 1
		}
	}
	if cull_bed != "" {
		if cull_bed, err = storage.Get(cull_bed, creds, false); err != nil {
			log.Printf("Error getting the cull file: %s\n", err.Error())
			return 1
		}
	}
	log.Printf("Time to get reference: %s\n", time.Since(ref_start))

	if use_metadata {
		mdata, err := extractMetadata(s1)
		if err != nil {
			log.Printf("Error extracting metadata from the file: %s\n", err.Error())
			return 1
		}
		sample_metadata = mdata
	} else {
		initializeMetadata(&sample_metadata, sample)
	}

	d1 = s1
	d2 = s2
	ramdisk = ""
	if use_ramdisk {
		ramdisk = "/dev/shm/"

		ram1 := ramdisk + common.GetLastSubstring(s1, "/")
		ram2 := ramdisk + common.GetLastSubstring(s2, "/")

		if err = common.CopyFile(s1, ram1); err != nil {
			log.Printf("Error loading sample (%s) into ramdisk (%s): %s\n", s1, ram1, err.Error())
			return 1
		}

		if err = common.CopyFile(s2, ram2); err != nil {
			log.Printf("Error loading sample (%s) into ramdisk (%s): %s\n", s2, ram2, err.Error())
			return 1
		}

		s1 = ram1
		s2 = ram2
	} else {
		//FIXME ramdisk is used as a global all over the place
		ramdisk = tmpDir
	}

	// Run alignment pipeline
	runPipe = losslessPipeline
	if lossy == true {
		runPipe = lossyPipeline
	}
	res_files, err = runPipe(s1, s2, ramdisk, !mark_repetitive, compression, tmpDir)
	if err != nil {
		log.Printf("Error in alignment pipeline: %s\n", err.Error())
		return 1
	}

	// at this point, the logs should be saved
	saveLogs = filepath.Join(op_path, sample+".log")

	// Picard -- PCR Duplicates
	if mark_repetitive == true {

		var runDup DuplicateOp

		runDup = picardMarkDuplicates
		if sambam == true {
			runDup = sambamMarkDuplicates
		}

		log.Printf("Marking PCR Duplicates ...")
		dup_start := time.Now()

		if dup_files, err := runDup(res_files[0], ramdisk, compression, tmpDir); err == nil {
			res_files = append(res_files, dup_files...)
		} else {
			log.Printf("Error marking PCR Duplicates: %s\n", err.Error())
			return 1
		}

		log.Printf("\t\t\t done: %s\n", time.Since(dup_start))
	}

	// Delete source files
	if remove_files == true {
		common.DeleteFile(d1)
		common.DeleteFile(d2)

		if use_ramdisk == true {
			common.DeleteFile(s1)
			common.DeleteFile(s2)
		}
	}

	// Write back sample and metadata
	if op_path != "" {

		log.Printf("Uploading samples to storage ...")
		upload_start := time.Now()
		for _, f := range res_files {
			log.Printf("\t\t\t%s --Copy-> %s", f, op_path)
			if err = storage.Put(op_path, f, creds, true); err != nil {
				log.Printf("Cannot save file (%s): %s\n", f, err.Error())
				return 1
			}
		}
		log.Printf("\t\t\t done: %s\n", time.Since(upload_start))
	} else {
		log.Printf("Files generated: %s\n", strings.Join(res_files, ", "))
	}

	log.Printf("Total runtime for %s: %s\n", sample, time.Since(ts_start))

	return 0
}

func main() {
	os.Exit(mainWithCode())
}
