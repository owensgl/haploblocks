package main

import "fmt"
import "time"
import "log"

import "genomics/common"

func ngmAlign(r1, r2, aligned string, isPaired, delOrig bool) (error) {

    // Multi-file: ngm -r <reference> -p --qry1 <r1> --qry2 <r2> -o <op> -t <nthreads>
    // Paired single: ngm -r <reference> -p -q <r1> -o <op> -t <nthreads>
    // Unpaired single: ngm -r <reference> -q <r1> -o <op> -t <nthreads>
    var ngm_args []string

    if r2 != "" && isPaired == false {
        return fmt.Errorf("Multiple files can't be unpaired reads")
    }

    log.Println("\tRunning ngm ...")
    ngm_args = []string{"-r", ref_genome, "-t", fmt.Sprintf("%d", nthreads), "-o", aligned}
    ngm_start := time.Now()

    if r2 != "" || isPaired == true {
        ngm_args = append(ngm_args, "-p")
    }

    if r2 != "" {
        ngm_args = append(ngm_args, "--qry1", r1, "--qry2", r2)
    } else {
        ngm_args = append(ngm_args, "-q", r1)
    }

    ngm_args = append(ngm_args, "--rg-id", sample_metadata.rg_id,
                                "--rg-pu", sample_metadata.rg_pu,
                                "--rg-sm", sample_metadata.real_name,
                                "--rg-lb", sample_metadata.real_name,
                                "--rg-pl", "ILLUMINA")

    if cmd_out, err := common.LogCommand(ngm_bin, ngm_args...).CombinedOutput(); err != nil {
        log.Printf("Error executing ngm: %s\n", string(cmd_out))
        return err
    }
    log.Printf("\t\t\t done: %s\n", time.Since(ngm_start))

    // Clean up input for paired files
    if delOrig == true {
        common.DeleteFile(r1)

        if r2 != "" {
            common.DeleteFile(r2)
        }
    }

    return nil
}
