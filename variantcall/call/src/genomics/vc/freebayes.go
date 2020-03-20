package main

import "fmt"
import "os"
import "genomics/common"

func freebayesVariants(bams []string, opfile, intervals string, hets float32, dummy_opts1 []string, dummy_opts2 []string, gvcf bool) (string, error) {

	fb_args := []string{"-f", ref_genome, "-v", opfile}

	if gvcf == true {
		fb_args = append(fb_args, "--gvcf")
	}

	if intervals != "" {
		fb_args = append(fb_args, "-t", intervals)
	}

	if len(bams) > 1 {
		lf := common.GetRandomString(8)

		f, err := os.Create(lf)
		if err != nil {
			return "", fmt.Errorf("Error opening output file: %s\n", err.Error())
		}
		defer common.DeleteFile(lf)

		for _, s := range bams {
			f.WriteString(s + "\n")
		}

		if err := f.Close(); err != nil {
			return "", fmt.Errorf("Error closing file: %s\n", err.Error())
		}

		fb_args = append(fb_args, "-L", lf)
	} else {
		fb_args = append(fb_args, bams[0])
	}

	cmd_out, err := common.LogCommand(vc_bin, fb_args...).CombinedOutput()
	if err != nil {
		newerr := fmt.Errorf("%s: freebayes: %s\n", err.Error(), string(cmd_out))
		return "", newerr
	}

	return opfile, nil
}
