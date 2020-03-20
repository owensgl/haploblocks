package main

import "fmt"
import "os/exec"
import "log"
import "os"
import "io"

func commonHaplotypeCaller(bams []string, opfile, intervals string, hets float32, java_options []string, vc_options []string,
	gvcf bool, gatkVersion string) (string, error) {

	var opflag string

	gatk_args := []string{}

	if java_options != nil && len(java_options) > 0 {
		gatk_args = append(gatk_args, java_options...)
	}

	gatk_args = append(gatk_args, "-jar", vc_bin)

	if gatkVersion == "4" || gatkVersion == "4beta" {
		gatk_args = append(gatk_args, "HaplotypeCaller")
		opflag = "-O"
	} else {
		gatk_args = append(gatk_args, "-T", "HaplotypeCaller")
		opflag = "-o"
	}
	gatk_args = append(gatk_args, "-R", ref_genome)

	for _, f := range bams {
		gatk_args = append(gatk_args, "-I", f)
	}
	if gvcf == true {
		if gatkVersion == "4" {
			gatk_args = append(gatk_args, "--emit-ref-confidence", "GVCF")
		} else {
			gatk_args = append(gatk_args, "--emitRefConfidence", "GVCF")
		}
	}

	if intervals != "" {
		gatk_args = append(gatk_args, "-L", intervals)
	}

	if gatkVersion == "4" {
		gatk_args = append(gatk_args, "--heterozygosity", fmt.Sprintf("%f", hets))

	} else {
		gatk_args = append(gatk_args, "-hets", fmt.Sprintf("%f", hets))
	}

	if vc_options != nil && len(vc_options) > 0 {
		gatk_args = append(gatk_args, vc_options...)
	}

	gatk_args = append(gatk_args, opflag, opfile)

	log.Printf("Invoking gatk command: java %s\n", gatk_args)
	log_file := opfile + ".log"

	log.Printf("Logs sent to: %s\n", log_file)

	cmd := exec.Command("java", gatk_args...)
	log_fd, err := os.OpenFile(log_file, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0755)
	if err != nil {
		log.Printf("Error opening log file (%s): %s\n", log_file, err.Error())
		return "", err
	}
	defer log_fd.Close()
	cmd.Stdout = log_fd
	cmd.Stderr = log_fd

	err = cmd.Run()
	if err != nil {
		// dump log to stderr
		_, errseek := log_fd.Seek(int64(0), 0)
		if errseek == nil {
			io.Copy(os.Stderr, log_fd)
		}
		newerr := fmt.Errorf("gatk haplotypecaller error: %s\n", err.Error())
		return "", newerr
	}

	return opfile, nil
}

func gatk3HaplotypeCaller(bams []string, opfile, intervals string, hets float32, java_options []string, extra_options []string, gvcf bool) (string, error) {
	return commonHaplotypeCaller(bams, opfile, intervals, hets, java_options, extra_options, gvcf, "3")
}

func gatk4BetaHaplotypeCaller(bams []string, opfile, intervals string, hets float32, java_options []string, extra_options []string, gvcf bool) (string, error) {
	return commonHaplotypeCaller(bams, opfile, intervals, hets, java_options, extra_options, gvcf, "4beta")
}

func gatk4HaplotypeCaller(bams []string, opfile, intervals string, hets float32, java_options []string, extra_options []string, gvcf bool) (string, error) {
	return commonHaplotypeCaller(bams, opfile, intervals, hets, java_options, extra_options, gvcf, "4")
}
