package bed

import "fmt"
import "strings"
import "os/exec"
import "log"

func MergeFbIntervals(bin string, intervals []string, ref_genome string, extra_opts []string, opfile string) (error) {

    bcf_args := []string{"concat"}

    bcf_args = append(bcf_args, intervals...)
    bcf_args = append(bcf_args, "-o", opfile)

    if extra_opts != nil && len(extra_opts) > 0 {
	    bcf_args = append(bcf_args, extra_opts...)
    }

    log.Printf("Invoking command: %s %s\n", bin, bcf_args)

    cmd_out, err := exec.Command(bin, bcf_args...).CombinedOutput()
    if err != nil {
        newerr := fmt.Errorf("%s: bcftools concat: %s\n", err.Error(), string(cmd_out))
        return newerr
    }

    return nil
}

func MergeGatk3Intervals(bin string, intervals []string, ref_genome string, extra_opts []string, opfile string) (error) {

    gatk_args := []string{"-cp", bin, "org.broadinstitute.gatk.tools.CatVariants",
                          "-R", ref_genome}

    for _, f := range(intervals) {
        gatk_args = append(gatk_args, "-V", f)
    }
    gatk_args = append(gatk_args, "-out", opfile)

    if extra_opts != nil && len(extra_opts) > 0 {
	    gatk_args = append(gatk_args, extra_opts...)
    }
    log.Printf("Invoking gatk command: java %s\n", gatk_args)

    cmd_out, err := exec.Command("java", gatk_args...).CombinedOutput()
    if err != nil {
        newerr := fmt.Errorf("%s: gatk3 CatVariants: %s\n", err.Error(), string(cmd_out))
        return newerr
    }

    return nil
}

func MergeGatk4Intervals(bin string, intervals []string, ref_genome string, extra_opts []string, opfile string) (error) {

    gatk_args := []string{"-jar", bin, "GatherVcfs"}

    for _, f := range(intervals) {
        gatk_args = append(gatk_args, "-I", f)
    }
    gatk_args = append(gatk_args, "-O", opfile)

    if extra_opts != nil && len(extra_opts) > 0 {
	    gatk_args = append(gatk_args, extra_opts...)
    }

    log.Printf("Invoking gatk command: java %s\n", gatk_args)

    cmd_out, err := exec.Command("java", gatk_args...).CombinedOutput()
    if err != nil {
        newerr := fmt.Errorf("%s: gatk4 GatherVcfs: %s\n", err.Error(), string(cmd_out))
        return newerr
    }

    return nil
}

func MergePlatypusIntervals(bin string, intervals []string, ref_genome string, extra_opts []string, opfile string) (error) {

    iplist := strings.Join(intervals, ",")
    plat_args := []string{bin, "mergeVariants",
                          fmt.Sprintf("--intervalFiles=%s", iplist),
                          fmt.Sprintf("--output=%s", opfile)}

    if extra_opts != nil && len(extra_opts) > 0 {
	    plat_args = append(plat_args, extra_opts...)
    }

    cmd_out, err := exec.Command("python", plat_args...).CombinedOutput()
    if err != nil {
        newerr := fmt.Errorf("%s: platypus concat: %s\n", err.Error(), string(cmd_out))
        return newerr
    }

    return nil
}
