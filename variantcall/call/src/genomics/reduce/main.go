package main

import "fmt"
import "flag"
import "strings"
import "log"
import "os"
import "bufio"
import "time"
import "math/rand"
import "path/filepath"

import "genomics/common"
import "genomics/bed"
import "genomics/storage"

var nthreads int
var ref_genome string
var vc_bin string

type MergeIntervalOp func(string, []string, string, []string, string) error

func getOutputName(subvcf string) string {

	ext := ".vcf"
	if strings.HasSuffix(subvcf, ".g.vcf") {
		ext = ".g.vcf"
	}

	idx := strings.Index(subvcf, ".bam")
	name := subvcf[:idx] + ext

	return name
}

func usage() {
	fmt.Println(
		` Usage: reduce -r REF -L LIST OPERATION CALLER [opts]*
	-r REF	Reference Genome
	-L LIS	Input List
	-u USR	 Username
	-p PWD	 Password
	-P NMS	 Project Namespace (GCP)
	-o OUT	Output Path

	-l	Log to file (default stdout)

	-vc PTH	Path to VC binary
	-m  PTH	Path to bcftools binary for merge

	OPERATION is exactly one of:

	-cat		(Concatenate different intervals into a single file).
	-genotype	(Perform a joint genotyping over individual samples)

	CALLER is exactly one of:

	-freebayes (Use freebayes)
	-platypus (Use platypus)
	-gatk3 (Use gatk3)
	-gatk4 (Use gatk4)
`)
}

func main() {

	const fb_path string = "/home/freebayes/freebayes"
	const plat_path string = "/home/platypus/Platypus.py"
	const gatk3_path string = "/home/gatk-3.8.0.jar"
	const gatk4_path string = "/gatk/gatk.jar"
	const bcf_path string = "/home/bcftools"
	const EINVAL = 0x16
	const ENOENT = 0x2
	const ENOTSUP = 0x5f

	var ref_path, ip_path string
	var ipfile string
	var vcffiles []string

	var cat, genotype bool

	var merge_path, merge_bin string
	var freebayes, platypus, gatk3, gatk4 bool

	var user, pass, project string
	var creds common.AuthenticationCreds
	var op_path, op_file string

	var to_log bool
	var log_file string
	var fd *os.File

	var runMerge MergeIntervalOp
	var ts_start, merge_start time.Time

	var err error

	ret := 0
	
	flag.StringVar(&ref_path, "r", "", "Reference Genome")
	flag.StringVar(&ip_path, "L", "", "Input List")

	flag.StringVar(&user, "u", "", "Username")
	flag.StringVar(&pass, "p", "", "Password")
	flag.StringVar(&project, "P", "", "Project Namespace (GCP)")

	flag.StringVar(&op_path, "o", "", "Output Path")
	flag.BoolVar(&to_log, "l", false, "Log to File")

	flag.BoolVar(&cat, "cat", false, "Concatenate")
	flag.BoolVar(&genotype, "genotype", false, "Genotype")

	flag.BoolVar(&freebayes, "freebayes", false, "Use freebayes")
	flag.BoolVar(&platypus, "platypus", false, "Use platypus")
	flag.BoolVar(&gatk3, "gatk3", false, "Use gatk3")
	flag.BoolVar(&gatk4, "gatk4", false, "Use gatk4")
	flag.StringVar(&merge_bin, "m", "", "Path to merge binary (freebayes only)")

	flag.Parse()
	
	rand.Seed(time.Now().UnixNano())
	if to_log == true {
		log_file = common.GetRandomString(8)
		if fd, err = os.Create(log_file); err != nil {
			panic(err)
		}
		log.SetOutput(fd)
	} else {
		log_file = ""
	}

	// Check all the files/inputs
	if ref_path == "" {
		log.Printf("Missing reference genome\n")
		usage()
		ret = -EINVAL
		goto out
	}

	if ip_path == "" {
		log.Printf("Missing input list\n")
		ret = -EINVAL
		goto out
	}

	if !(cat || genotype) {
		log.Printf("No operation specified\n")
		usage()
		ret = -EINVAL
		goto out
	}

	if cat && genotype {
		log.Printf("Cannot concatenate and genotype\n")
		usage()
		ret = -EINVAL
		goto out
	}

	if !(freebayes || platypus || gatk3 || gatk4) {
		log.Printf("No variant caller specified\n")
		usage()
		ret = -EINVAL
		goto out
	}

	if (freebayes || platypus || gatk3) &&
		(platypus || gatk3 || gatk4) &&
		(gatk3 || gatk4 || freebayes) &&
		(gatk4 || freebayes || platypus) {
		log.Printf("Multiple variant callers specified\n")
		usage()
		ret = -EINVAL
		goto out
	}

	if (gatk3 && merge_bin != "") || (gatk4 && merge_bin != "") ||
		(platypus && merge_bin != "") {
		log.Printf("Merge option only valid with -freebayes\n")
		usage()
		ret = -EINVAL
		goto out
	}

	if genotype {
		log.Printf("Genotyping currently unsupported\n")
		ret = -ENOTSUP
		goto out
	}

	if (user != "" && pass == "") || (user == "" && pass != "") {
		log.Printf("Incomplete credentials\n")
		usage()
		goto out
	}
	creds = common.AuthenticationCreds{Username: user, Password: pass, Project: project}

	// Acquire files
	ts_start = time.Now()
	if ref_genome, err = storage.Get(ref_path, creds, false); err != nil {
		log.Printf("Error getting the reference: %s\n", err.Error())
		ret = -ENOENT
		goto out
	}

	if ipfile, err = storage.Get(ip_path, creds, false); err != nil {
		log.Printf("Error getting the input list: %s\n", err.Error())
		ret = -ENOENT
		goto out
	}

	fd, err = os.Open(ipfile)
	if err == nil {
		defer fd.Close()

		scanner := bufio.NewScanner(fd)
		for scanner.Scan() {

			line := scanner.Text()
			f, e := storage.Get(line, creds, false)
			if e != nil {
				log.Printf("Error getting %s: %s\n", line, e.Error())
				goto out
			}

			vcffiles = append(vcffiles, f)
		}

	} else {
		log.Printf("Error opening file (%s): %s\n", ipfile, err.Error())
		ret = -ENOENT
		goto out
	}
	log.Printf("Time to get files: %s\n", time.Since(ts_start))

	// Run selected variant caller
	if freebayes == true {
		merge_path = bcf_path
		runMerge = bed.MergeFbIntervals
	} else if platypus == true {
		merge_path = plat_path
		merge_bin = vc_bin
		runMerge = bed.MergePlatypusIntervals
	} else if gatk3 == true {
		merge_path = gatk3_path
		merge_bin = vc_bin
		runMerge = bed.MergeGatk3Intervals
	} else if gatk4 == true {
		merge_path = gatk4_path
		merge_bin = vc_bin
		runMerge = bed.MergeGatk4Intervals
	}

	if merge_bin == "" {
		merge_bin = merge_path
	}

	log.Printf("Merging files ...\n")
	op_file = getOutputName(vcffiles[0])
	merge_start = time.Now()
	if err = runMerge(merge_bin, vcffiles, ref_genome, nil, op_file); err != nil {
		log.Printf("Error merging the intervals: %s\n", err.Error())
		ret = 1
		goto out
	}
	log.Printf("\t\t\t done: %s\n", time.Since(merge_start))

	for _, f := range vcffiles {
		common.DeleteFile(f)
		if common.FileExists(f + ".idx") {
			common.DeleteFile(f + ".idx")
		}
	}

	// Write back sample and metadata
	if op_path != "" {

		log.Println("Uploading to storage ...")
		upload_start := time.Now()
		if err = storage.Put(op_path, op_file, creds, true); err != nil {
			log.Printf("Cannot save file (%s): %s\n", op_file, err.Error())
			ret = 1
			goto out
		}

		idx_file := op_file + ".idx"
		if common.FileExists(idx_file) {
			if err = storage.Put(op_path, idx_file, creds, true); err != nil {
				log.Printf("Cannot save index file (%s): %s\n", idx_file, err.Error())
				ret = 1
				goto out
			}
		}

		log.Printf("\t\t\t done: %s\n", time.Since(upload_start))
	} else {
		log.Printf("Files generated: %s\n", op_file)
	}

	log.Printf("Total runtime for variant calling: %s\n", time.Since(ts_start))
	ret = 0
out:
	if log_file != "" {
		fd.Sync()
		fd.Close()

		up_log := filepath.Base(ipfile) + ".log"
		os.Rename(log_file, up_log)

		if op_path != "" {
			storage.Put(op_path, up_log, creds, true)
		}
	}
	os.Exit(ret)
}
