package main

import "fmt"
import "strings"
import "genomics/common"

func platypusVariants(bams []string, opfile, intervals string, hets float32, dummy_opts1 []string, dummy_opts2 []string, gvcf bool) (string, error) {

	iplist := strings.Join(bams, ",")
	plat_args := []string{vc_bin, "callVariants",
		fmt.Sprintf("--nCPU=%d", 1),
		fmt.Sprintf("--refFile=%s", ref_genome),
		fmt.Sprintf("--bamFiles=%s", iplist),
		fmt.Sprintf("--output=%s", opfile)}

	if gvcf == true {
		plat_args = append(plat_args, "--outputRefCalls=1", "--refCallBlockSize=128")
	}

	if intervals != "" {
		plat_args = append(plat_args, fmt.Sprintf("--regions=%s", intervals))
	}

	cmd_out, err := common.LogCommand("python", plat_args...).CombinedOutput()
	if err != nil {
		newerr := fmt.Errorf("%s: platypus: %s\n", err.Error(), string(cmd_out))
		return "", newerr
	}

	return opfile, nil
}
