package main

import "fmt"
import "log"
import "os"
import "path/filepath"
import "genomics/common"

func samToBam(sam, bam string, delOrig bool) error {

	// samtools view -Su <sam> -o <bam>
	convert_args := []string{"view", "-Su", sam, "-o", bam}

	cmd_out, err := common.LogCommand("samtools", convert_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in sam to bam conversion: %s\n", string(cmd_out))
		return err
	}

	if delOrig == true {
		common.DeleteFile(sam)
	}

	return nil
}

func samSort(file string, comp int, tmpDir string) error {

	// samtools sort -f -l 0 -m <heap> -@ <nthreads> <ip> <op>
	// Pick a heap of 8G
	tmpfile := filepath.Dir(file) + "/" + common.GetRandomString(16)
	os.Rename(file, tmpfile)
	defer common.DeleteFile(tmpfile)

	sort_args := []string{"sort", "-f", "-l", fmt.Sprintf("%d", comp),
		"-m", "8000000000",
		"-@", fmt.Sprintf("%d", nthreads),
		tmpfile, file}

	cmd_out, err := common.LogCommand("samtools", sort_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in samtools sort: %s\n", string(cmd_out))
	}

	return err
}

func samIndex(file string) (string, error) {
	index_args := []string{"index", file}

	cmd_out, err := common.LogCommand("samtools", index_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in samtools index: %s\n", string(cmd_out))
		return "", err
	}

	return file + ".bai", nil
}
