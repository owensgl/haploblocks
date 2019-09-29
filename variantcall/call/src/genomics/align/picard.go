package main

import "fmt"
import "log"
import "os"
import "path/filepath"

import "genomics/common"

func picardFastqToSam(fastq1, fastq2, sam string, comp int, sort bool,
	real_name, read_group, platform_unit string, delOrig bool) error {

	default_sort := "queryname"
	if sort == false {
		default_sort = "unsorted"
	}

	picard_args := []string{"-jar", picard_bin, "FastqToSam",
		fmt.Sprintf("FASTQ=%s", fastq1),
		fmt.Sprintf("FASTQ2=%s", fastq2),
		fmt.Sprintf("OUTPUT=%s", sam),
		fmt.Sprintf("SORT_ORDER=%s", default_sort),
		fmt.Sprintf("COMPRESSION_LEVEL=%d", comp),
		fmt.Sprintf("READ_GROUP_NAME=%s", read_group),
		fmt.Sprintf("SAMPLE_NAME=%s", real_name),
		fmt.Sprintf("LIBRARY_NAME=%s", real_name),
		fmt.Sprintf("PLATFORM_UNIT=%s", platform_unit),
		"PLATFORM=illumina",
		"SEQUENCING_CENTER=GenomeQuebec"}

	cmd_out, err := common.LogCommand("java", picard_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in picard fastqtosam: %s\n", string(cmd_out))
		return err
	}

	if delOrig == true {
		common.DeleteFile(fastq1)
		common.DeleteFile(fastq2)
	}

	return nil
}

func picardSamToFastq(sam, fastq string, delOrig bool) error {

	picard_args := []string{"-jar", picard_bin, "SamToFastq",
		fmt.Sprintf("INPUT=%s", sam),
		fmt.Sprintf("FASTQ=%s", fastq),
		"CLIPPING_ATTRIBUTE=XT",
		"CLIPPING_ACTION=2",
		"INTERLEAVE=true",
		"NON_PF=true"}

	cmd_out, err := common.LogCommand("java", picard_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in picard samtofastq: %s\n", string(cmd_out))
		return err
	}

	if delOrig == true {
		common.DeleteFile(sam)
	}

	return nil
}

func picardSortEx(ip, op string, comp int) error {

	picard_args := []string{"-jar", picard_bin, "SortSam",
		fmt.Sprintf("I=%s", ip),
		fmt.Sprintf("O=%s", op),
		fmt.Sprintf("COMPRESSION_LEVEL=%d", comp),
		"SORT_ORDER=queryname"}

	cmd_out, err := common.LogCommand("java", picard_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in picard sort: %s\n", string(cmd_out))
	}

	return err
}

func picardSort(file string, comp int) error {

	tmpfile := filepath.Dir(file) + "/" + common.GetRandomString(16)
	os.Rename(file, tmpfile)
	defer common.DeleteFile(tmpfile)

	return picardSortEx(tmpfile, file, comp)
}

func picardMarkIlluminaAdapters(file string, output_dir string, comp int) (string, string, error) {

	tmpfile := filepath.Dir(file) + "/" + common.GetRandomString(16)
	metfile := filepath.Join(output_dir, sample) + ".illuminametrics.txt"

	os.Rename(file, tmpfile)
	defer common.DeleteFile(tmpfile)

	picard_args := []string{"-jar", picard_bin, "MarkIlluminaAdapters",
		fmt.Sprintf("I=%s", tmpfile),
		fmt.Sprintf("O=%s", file),
		fmt.Sprintf("M=%s", metfile),
		fmt.Sprintf("COMPRESSION_LEVEL=%d", comp)}

	cmd_out, err := common.LogCommand("java", picard_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in picard markilluminaadapters: %s\n", string(cmd_out))
		return "", "", err
	}

	return file, metfile, nil
}

func picardMergeBamAlignment(bam, ubam, merged, reference string,
	comp int, create_index, create_hash, delOrig bool) ([]string, error) {

	md5file := merged + ".md5"
	indexfile := merged[:len(merged)-1] + "i"

	picard_args := []string{"-jar", picard_bin, "MergeBamAlignment",
		fmt.Sprintf("ALIGNED=%s", bam),
		fmt.Sprintf("UNMAPPED=%s", ubam),
		fmt.Sprintf("OUTPUT=%s", merged),
		fmt.Sprintf("R=%s", reference),
		fmt.Sprintf("COMPRESSION_LEVEL=%d", comp),
		fmt.Sprintf("CREATE_INDEX=%t", create_index),
		fmt.Sprintf("CREATE_MD5_FILE=%t", create_hash),
		"CLIP_ADAPTERS=false",
		"CLIP_OVERLAPPING_READS=true",
		"MAX_INSERTIONS_OR_DELETIONS=-1"}

	cmd_out, err := common.LogCommand("java", picard_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in picard mergebamalignment: %s\n", string(cmd_out))
		return nil, err
	}

	if delOrig == true {
		common.DeleteFile(ubam)
		common.DeleteFile(bam)
	}

	res_files := []string{merged}
	if create_index == true {
		res_files = append(res_files, indexfile)
	}
	if create_hash == true {
		res_files = append(res_files, md5file)
	}

	return res_files, nil
}

func picardMarkDuplicates(file string, output_dir string, comp int, tmpDir string) ([]string, error) {

	tmpfile := filepath.Join(tmpDir, "tmp-markdup-input-"+common.GetRandomString(16))
	metfile := filepath.Join(output_dir, sample) + ".dupmetrics.txt"
	md5file := file + ".md5"
	indexfile := file[:len(file)-1] + "i"

	os.Rename(file, tmpfile)
	defer common.DeleteFile(tmpfile)

	picard_args := []string{"-jar", picard_bin, "MarkDuplicates",
		fmt.Sprintf("I=%s", tmpfile),
		fmt.Sprintf("O=%s", file),
		fmt.Sprintf("M=%s", metfile),
		fmt.Sprintf("COMPRESSION_LEVEL=%d", comp),
		fmt.Sprintf("TMP_DIR=%s", tmpDir),
		"VALIDATION_STRINGENCY=LENIENT",
		"CREATE_INDEX=TRUE",
		"CREATE_MD5_FILE=TRUE"}

	cmd_out, err := common.LogCommand("java", picard_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in picard markdups: %s\n", string(cmd_out))
		return nil, err
	}

	return []string{metfile, indexfile, md5file}, nil
}
