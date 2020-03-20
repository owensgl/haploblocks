package main

import "fmt"
import "log"
import "path/filepath"
import "genomics/common"

func trimmomatic(s1, s2, output_path string, trim_prefix string, delOrig bool) (string,
	string, string, string, error) {

	// Drop poor quality reads
	// java -jar trimmomatic.jar PE -threads <count> -phred33 <fwd> <rev>
	// <fwd_pair> <fwd_unpair> <rev_pair> <rev_unpair>
	// ILLUMINACLIP:<sample_prefix>.fa:2:30:10:8:T SLIDINGWINDOW:4:15 MINLEN:36

	s1_pair := filepath.Join(output_path, sample+"_R1.paired.fastq")
	s1_unpair := filepath.Join(output_path, sample+"_R1.unpaired.fastq")
	s2_pair := filepath.Join(output_path, sample+"_R2.paired.fastq")
	s2_unpair := filepath.Join(output_path, sample+"_R2.unpaired.fastq")

	trim_args := []string{"-jar", trimmomatic_bin, "PE",
		"-threads", fmt.Sprintf("%d", nthreads), "-phred33",
		s1, s2, s1_pair, s1_unpair, s2_pair, s2_unpair,
		fmt.Sprintf("ILLUMINACLIP:%s:2:30:10:8:T", trim_prefix),
		"SLIDINGWINDOW:4:15", "MINLEN:36"}

	cmd_out, err := common.LogCommand("java", trim_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error trimming %s and %s: %s\n", s1, s2, string(cmd_out))
		return "", "", "", "", err
	}

	// Clean heck to see if all files exist and clean up
	if delOrig == true {
		common.DeleteFile(s1)
		common.DeleteFile(s2)
	}

	return s1_pair, s2_pair, s1_unpair, s2_unpair, nil
}

func cullNonrepetitiveRegions(file, cull_bed, culled_file string) error {

	// Remove non-repetitive regions
	log.Printf("Remove non-repetitive regions ...")
	cull_args := []string{"writeRegion", "--bed", cull_bed, "--in", file,
		"--out", culled_file}

	cmd_out, err := common.LogCommand(bamutils_bin, cull_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in bamutils writeregion: %s\n", string(cmd_out))
	}

	return err
}
