package main

import "fmt"
import "strings"
import "flag"
import "log"
import "os"
import "io/ioutil"
import "bufio"
import "time"
import "runtime"
import "math/rand"
import "path/filepath"
import "sync"

import "genomics/common"
import "genomics/bed"
import "genomics/storage"

var nthreads int
var plat_threads int
var ref_genome string
var vc_bin string
var log_tmp string
var log_fd *os.File
var creds common.AuthenticationCreds

/* bamfiles, opfile, req.intervalFile, 0.01, java_options, vc_options, gvcf */
type CallVariantOp func([]string, string, string, float32, []string, []string, bool) (string, error)

/* merge_bin, resfiles, ref_genome, merge_options, op_file */
type MergeIntervalOp func(string, []string, string, []string, string) error

func getOutputName(bam, tmpDir string, intervals string, gvcf bool) string {

	const marker string = "intervals"

	opfile := filepath.Join(tmpDir, filepath.Base(bam))

	/// TODO: Need to encode arbitrary bed names into the result too
	if strings.Contains(intervals, marker) {
		i := strings.LastIndex(intervals, marker) + len(marker)
		opfile += intervals[i:(len(intervals) - 4)]
	}

	if gvcf == true {
		opfile += ".g.vcf"
	} else {
		opfile += ".vcf"
	}

	return opfile
}

func usage() {
	fmt.Println("Usage:")
	fmt.Println("\t -n \t Thread Count")
	fmt.Println("\t -r \t Reference Genome")
	fmt.Println("\t -i \t Input File")
	fmt.Println("\t -L \t Input List")
	fmt.Println("\t -b PATH\t Intervals to Call (BED format)")
	fmt.Println("\t -g \t Generate GVCF (default false)")
	fmt.Println("\t -u USER\t Username")
	fmt.Println("\t -p PASS\t Password")
	fmt.Println("\t -P NAME\t Project Name (GCP)")
	fmt.Println("\t -o PATH\t Output Path")
	fmt.Println("\t -l PATH\t Log to File (default stdout)")
	fmt.Println("\t -w DIR\t Work directory (for temporary files and such. defaults to '.')")
	fmt.Println("\t    \t Place on the same filesystem as output path for efficiency.")
	fmt.Println("\t -nsegments NUM\t Number of segments in which the ")
	fmt.Println("\t -minbp NUM\t Only consider contigs from BED which are >= NUM bp")
	fmt.Println("\t -----------------------------------")
	fmt.Println("\t -freebayes (Use freebayes)")
	fmt.Println("\t -platypus (Use platypus)")
	fmt.Println("\t -gatk3 (Use gatk3)")
	fmt.Println("\t -gatk4beta (Use gatk4 beta)")
	fmt.Println("\t -gatk4 (Use gatk4)")
	fmt.Println("\t -vc \t Path to VC binary")
	fmt.Println("\t -m \t Path to bcftools binary for merge")

	os.Exit(0)
}

func closeLogs(savePath *string, uploadPath *string) {
	if log_tmp != "" {
		if log_fd != nil {
			log_fd.Sync()
			log_fd.Close()
		}

		if *savePath != "" {
			os.Rename(log_tmp, *savePath)
			if *uploadPath != "" {
				storage.Put(*uploadPath, *savePath, creds, true)
			}
		} else {
			os.Remove(log_tmp)
		}
	}
}

type VCWorkRequest struct {
	intervalNo   int
	intervalFile string
	deleteBed    bool
}

type VCWorkResult struct {
	intervalNo int
	err        error
	resFile    string
}

// can't believe golang needs this
func bToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

func mainWithCode() int {

	const fb_path string = "/home/freebayes/freebayes"
	const plat_path string = "/home/platypus/Platypus.py"
	const gatk3_path string = "/home/gatk-3.8.0.jar"
	const gatk4_path string = "/gatk/gatk.jar"
	const bcf_path string = "/home/bcftools"

	var ref_path, bam_path, bam_list, bed_path string

	var java_options_s string
	var vc_options_s string
	var merge_options_s string
	java_options  := []string{}
	vc_options    := []string{}
	merge_options := []string{}

	var bamfile, bedfile string
	var bamfiles []string
	var intervals []string
	var resfiles []string

	var gvcf bool
	var delete_beds bool

	var vc_path, merge_path, merge_bin string
	var freebayes, platypus, gatk3, gatk4, gatk4beta bool

	var user, pass, project string
	var op_path, op_file string

	var to_log bool
	var fd *os.File

	var runVC CallVariantOp
	var runMerge MergeIntervalOp

	var ts_start, vc_start time.Time

	var err error
	var barrier sync.WaitGroup
	var workDir string
	var credsFile string
	var nsegments int
	var minbp uint64

	flag.IntVar(&nthreads, "n", runtime.NumCPU(), "Thread Count")
	flag.IntVar(&nsegments, "nsegments", 0, "Number of segments (defaults to 1 per thread)")
	flag.Uint64Var(&minbp, "minbp", uint64(0), "Only include contigs in BED >= N bp. (defaults to 0, i.e. all of BED)")
	flag.StringVar(&ref_path, "r", "", "Reference Genome")
	flag.StringVar(&bam_path, "i", "", "Input File")
	flag.StringVar(&bam_list, "L", "", "Input List")
	flag.StringVar(&bed_path, "b", "", "Intervals to Call (BED format)")
	flag.BoolVar(&gvcf, "g", false, "Generate GVCF (default false)")

	flag.StringVar(&user, "u", "", "Username")
	flag.StringVar(&pass, "p", "", "Password")
	flag.StringVar(&project, "P", "", "Project Namespace (GCP)")
	flag.StringVar(&java_options_s,  "javaoptions",   "", "Additional options to pass to java command line. Pass multiple separated by spaces.")
	flag.StringVar(&vc_options_s,    "vcoptions",    "", "Additional options to pass to the vc tool.")
	flag.StringVar(&merge_options_s, "mergeoptions", "", "Additional options to pass to the merge tool.")

	flag.StringVar(&op_path, "o", "", "Output Path")
	flag.BoolVar(&to_log, "l", false, "Log to File")
	flag.StringVar(&workDir, "w", ".", "Work directory")
	flag.StringVar(&vc_bin, "vc", "", "Path to VC binary")
	flag.StringVar(&merge_bin, "m", "", "Path to merge binary (freebayes only)")
	flag.StringVar(&credsFile, "creds", "", "Path to credentials file")
	flag.BoolVar(&freebayes, "freebayes", false, "Use freebayes")
	flag.BoolVar(&platypus, "platypus", false, "Use platypus")
	flag.BoolVar(&gatk3, "gatk3", false, "Use gatk3")
	flag.BoolVar(&gatk4, "gatk4", false, "Use gatk4")
	flag.BoolVar(&gatk4beta, "gatk4beta", false, "Use gatk4 beta")

	flag.Parse()

	rand.Seed(time.Now().UnixNano())
	if to_log == true {
		log_tmp = common.GetRandomString(8)
		if fd, err = os.Create(log_tmp); err != nil {
			panic(err)
		}
		log.SetOutput(fd)
	} else {
		log_tmp = ""
	}

	local_log := "vc.log"
	defer closeLogs(&local_log, &op_path)

	// Check all the files/inputs
	if ref_path == "" {
		log.Printf("Missing reference genome\n")
		usage()
		return 1
	}

	if bam_path == "" && bam_list == "" {
		log.Printf("Missing input file\n")
		usage()
		return 1
	}
	if bam_path != "" {
		local_log = filepath.Base(bam_path) + ".log"
		if op_path == "" {
			op_path = filepath.Dir(bam_path)
		}
	} else {
		local_log = filepath.Base(bam_list) + ".log"
		if op_path == "" {
			op_path = filepath.Dir(bam_list)
		}
	}

	if bam_path != "" && bam_list != "" {
		log.Printf("Can't specify both input file and list\n")
		usage()
		return 1
	}

	if (user != "" && pass == "") || (user == "" && pass != "") {
		log.Printf("Incomplete credentials\n")
		usage()
		return 1
	}

	creds, err = common.ParseCreds(user, pass, project, credsFile)
	if err != nil {
		log.Printf("can't parse credentials: %s", err)
		return 1
	}

	var numCallers = bToInt(freebayes) + bToInt(platypus) + bToInt(gatk4) + bToInt(gatk3) + bToInt(gatk4beta)

	if numCallers == 0 {
		log.Printf("No variant caller specified\n")
		usage()
		return 1
	}

	if numCallers > 1 {
		log.Printf("Multiple variant callers specified\n")
		usage()
		return 1
	}

	if (gatk3 && merge_bin != "") || (gatk4 && merge_bin != "") ||
		(platypus && merge_bin != "") || (gatk4beta && merge_bin != "") {
		log.Printf("Merge option only valid with -freebayes\n")
		usage()
		return 1
	}

	// Setup working directory to store all files
	var info os.FileInfo
	info, err = os.Stat(workDir)
	if err != nil || !info.IsDir() {
		log.Printf("invalid working directory: %s", workDir)
		return 1
	}
	tmpDir, err := ioutil.TempDir(workDir, ".vc-")
	if err != nil {
		log.Printf("could not create temp directory inside %s: %s", workDir, err.Error())
		return 1
	}
	storage.LocalDir = tmpDir
	defer os.RemoveAll(tmpDir)

	// Acquire files
	ts_start = time.Now()
	if ref_genome, err = storage.Get(ref_path, creds, false); err != nil {
		log.Printf("Error getting the reference: %s\n", err.Error())
		return 1
	}

	if bam_path != "" {
		if bamfile, err = storage.Get(bam_path, creds, false); err != nil {
			log.Printf("Error getting the input file: %s\n", err.Error())
			return 1
		}
		bamfiles = make([]string, 1, 1)
		bamfiles[0] = bamfile
	}

	if bam_list != "" {
		if bamfile, err = storage.Get(bam_list, creds, false); err != nil {
			log.Printf("Error getting the input list: %s\n", err.Error())
			return 1
		}

		fd, e := os.Open(bamfile)
		if e != nil {
			log.Printf("Error opening file (%s): %s\n", err.Error())
			return 1
		}
		defer fd.Close()

		scanner := bufio.NewScanner(fd)
		for scanner.Scan() {

			line := scanner.Text()
			f, e := storage.Get(line, creds, false)
			if e != nil {
				log.Printf("Error getting %s: %s\n", line, err.Error())
				return 1
			}

			bamfiles = append(bamfiles, f)
		}
	}

	if bed_path != "" {
		if bedfile, err = storage.Get(bed_path, creds, false); err != nil {
			log.Printf("Error getting the interval file: %s\n", err.Error())
			return 1
		}
	} else {
		log.Printf("Missing interval file\n")
		return 1
	}
	log.Printf("Time to get files: %s\n", time.Since(ts_start))

	// Run selected variant caller
	if freebayes == true {
		vc_path = fb_path
		merge_path = bcf_path
		runVC = freebayesVariants
		runMerge = bed.MergeFbIntervals
	} else if platypus == true {
		vc_path = plat_path
		merge_path = plat_path
		merge_bin = vc_bin
		runVC = platypusVariants
		runMerge = bed.MergePlatypusIntervals
	} else if gatk3 == true {
		vc_path = gatk3_path
		merge_path = gatk3_path
		merge_bin = vc_bin
		runVC = gatk3HaplotypeCaller
		runMerge = bed.MergeGatk3Intervals
	} else if gatk4beta == true {
		vc_path = gatk4_path
		merge_path = gatk4_path
		merge_bin = vc_bin
		runVC = gatk4BetaHaplotypeCaller
		runMerge = bed.MergeGatk4Intervals

		/* XXX: GATK4 handles parallelism internally using Spark runner.
		 * At the moment, the documentation seems pretty unclear about
		 * calling Spark vs non-spark versions of tools. Relying wholly
		 * on the runner means CPUs are still idle around 50%;
		 * giving it the full core count actually slows stuff done.
		 * Picking a magic number of /4. Can tweak later.
		 */
		nthreads = (nthreads + 3) / 4
	} else if gatk4 == true {
		vc_path = gatk4_path
		merge_path = gatk4_path
		merge_bin = vc_bin
		runVC = gatk4HaplotypeCaller
		runMerge = bed.MergeGatk4Intervals

		/* XXX: GATK4 handles parallelism internally using Spark runner.
		 * At the moment, the documentation seems pretty unclear about
		 * calling Spark vs non-spark versions of tools. Relying wholly
		 * on the runner means CPUs are still idle around 50%;
		 * giving it the full core count actually slows stuff done.
		 * Picking a magic number of /4. Can tweak later.
		 */
		nthreads = (nthreads + 3) / 4
	}

	if vc_bin == "" {
		vc_bin = vc_path
	}

	if merge_bin == "" {
		merge_bin = merge_path
	}

	if nsegments == 0 {
		nsegments = nthreads
	}

	if nsegments < nthreads {
		nthreads = nsegments
	}

	if nsegments > 1 {
		log.Printf("Splitting bed into %d segments\n", nsegments)
		pid := fmt.Sprintf("%d", os.Getpid())
		dstPrefix := filepath.Join(tmpDir, pid)
		if intervals, err = bed.SplitIntervals(bedfile, dstPrefix, nsegments, minbp); err != nil {
			log.Printf("Error in splitting the bed file: %s\n", err.Error())
			return 1
		}
		for _, interval := range intervals {
			log.Printf("  => %s", interval)
		}
		delete_beds = false
	} else {
		intervals = make([]string, 1, 1)
		intervals[0] = bedfile
		delete_beds = false
	}


	if java_options_s != "" {
		xopts := strings.Split(java_options_s, " ")
		java_options = append(java_options, xopts...)
	}
	if vc_options_s != "" {
		xopts := strings.Split(vc_options_s, " ")
		vc_options = append(vc_options, xopts...)
	}
	if merge_options_s != "" {
		xopts := strings.Split(merge_options_s, " ")
		merge_options = append(merge_options, xopts...)
	}

	log.Printf("Variant calling with %d threads, %d segments...\n", nthreads, nsegments)

	var workQueue chan VCWorkRequest = make(chan VCWorkRequest, nsegments)
	var resultQueue chan VCWorkResult = make(chan VCWorkResult, nsegments)

	resfiles = make([]string, nsegments, nsegments)

	barrier.Add(nthreads)
	vc_start = time.Now()

	// start workers
	for i := 0; i < nthreads; i++ {
		go func(wq chan VCWorkRequest, rq chan VCWorkResult, workerId int) {
			defer barrier.Done()

			for {
				req, ok := <-wq
				if !ok {
					log.Printf("\t[worker %d] Tearing down...", workerId)
					break
				}

				segment_start := time.Now()
				opfile := getOutputName(bamfile, tmpDir, req.intervalFile, gvcf)
				resfile, err := runVC(bamfiles, opfile, req.intervalFile, 0.01, java_options, vc_options, gvcf)
				result := VCWorkResult{
					intervalNo: req.intervalNo,
					err:        err,
					resFile:    resfile}

				if err != nil {
					log.Printf("\t[worker %d] Error in alignment pipeline: %s\n", workerId, err.Error())

				} else {
					log.Printf("\t[worker %d] Segment %d complete: %s\n", workerId, req.intervalNo, time.Since(segment_start))
				}

				if req.deleteBed {
					common.DeleteFile(req.intervalFile)
					if gatk3 && err != nil {
						common.DeleteFile(resfile + ".idx")
					}
				}
				rq <- result
			}
		}(workQueue, resultQueue, i)
	}

	// submit some work ahead of time
	for i := 0; i < nthreads; i++ {
		workQueue <- VCWorkRequest{
			intervalNo:   i,
			intervalFile: intervals[i],
			deleteBed:    delete_beds}
	}

	nextSegment := nthreads
	segmentDone := 0
	problems := 0
	for {
		res := <-resultQueue
		if res.err != nil {
			log.Printf("Error occurred in one worker. Tearing down other workers...")
			problems++
			break
		}

		resfiles[res.intervalNo] = res.resFile

		segmentDone++
		if segmentDone >= nsegments {
			log.Printf("All segments done. tearing down.")
			break
		} else {
			log.Printf("Completed %d/%d segments (%d%%): %s", segmentDone, nsegments, segmentDone*100.0/nsegments, time.Since(vc_start))
		}

		if nextSegment < nsegments {
			workQueue <- VCWorkRequest{
				intervalNo:   nextSegment,
				intervalFile: intervals[nextSegment],
				deleteBed:    delete_beds}
			nextSegment++
		}
	}
	close(workQueue)
	barrier.Wait()

	log.Printf("Runtime for variant calling on all segments: %s\n", time.Since(vc_start))
	if problems > 0 {
		log.Printf("Problems occurred. Cannot merge.")
		return 1
	}

	if nsegments > 1 {

		log.Printf("Merging interval files ...\n")
		op_file = getOutputName(bamfile, tmpDir, bedfile, gvcf)

		// Check to see interval vcfs exist
		for _, f := range resfiles {
			if common.FileExists(f) == false {
				log.Printf("Missing interval file %s, exiting\n", f)
				return 1
			}
		}

		merge_start := time.Now()
		if err = runMerge(merge_bin, resfiles, ref_genome, merge_options, op_file); err != nil {
			log.Printf("Error merging the intervals: %s\n", err.Error())
			return 1
		}
		log.Printf("\t\t\t done: %s\n", time.Since(merge_start))

		if !platypus {
			for _, f := range resfiles {
				common.DeleteFile(f)
			}
		}
	} else {
		op_file = resfiles[0]
	}

	// Write back sample and index
	if op_path != "" {

		genfiles := []string{op_file}
		idx_file := op_file + ".idx"
		if common.FileExists(idx_file) {
			genfiles = append(genfiles, idx_file)
		}

		log.Println("Uploading to storage ...")
		upload_start := time.Now()

		for _, f := range genfiles {
			log.Printf("\t\t\t uploading: %s => %s\n", f, op_path)
			if err = storage.Put(op_path, f, creds, true); err != nil {
				log.Printf("Cannot save file (%s): %s\n", f, err.Error())
				return 1
			}
			common.DeleteFile(f)
		}

		log.Printf("\t\t\t done: %s\n", time.Since(upload_start))
	} else {
		log.Printf("Files generated: %s\n", op_file)
	}

	log.Printf("Total runtime for variant calling: %s\n", time.Since(ts_start))
	return 0
}

func main() {
	os.Exit(mainWithCode())
}
