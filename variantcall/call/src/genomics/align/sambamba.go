package main

import "fmt"
import "log"
import "os"
import "path/filepath"

import "genomics/common"

func sambamSortEx(file string, comp int, picard bool, tmpDir string) error {

	// sambamba sort -l 0 -t <nthreads> -n -o <op> <ip>
	tmpfile := filepath.Dir(file) + "/" + common.GetRandomString(16)
	os.Rename(file, tmpfile)
	defer common.DeleteFile(tmpfile)

	sort_args := []string{"sort", "-l", fmt.Sprintf("%d", comp),
		"-t", fmt.Sprintf("%d", nthreads)}

	if picard == true {
		sort_args = append(sort_args, "-n")
	}

	sort_args = append(sort_args, "--tmpdir="+tmpDir, "-o", file, tmpfile)

	cmd_out, err := common.LogCommand(sambamba_bin, sort_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in sambamba sort: %s\n", string(cmd_out))
	}

	return err
}

func sambamSort(file string, comp int, tmpDir string) error {
	return sambamSortEx(file, comp, false, tmpDir)
}

func sambamSortPicard(file string, comp int, tmpDir string) error {
	return sambamSortEx(file, comp, true, tmpDir)
}

func sambamIndex(file string) (string, error) {

	indexfile := file + ".bai"
	index_args := []string{"index", "-t", fmt.Sprintf("%d", nthreads),
		file, indexfile}

	cmd_out, err := common.LogCommand(sambamba_bin, index_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in sambamba index: %s\n", string(cmd_out))
		return "", err
	}

	return indexfile, nil
}

func sambamMarkDuplicates(file string, output_dir string, comp int, tmpDir string) ([]string, error) {

	tmpfile := filepath.Join(tmpDir, "tmp-markdup-input-"+common.GetRandomString(16))

	os.Rename(file, tmpfile)
	defer common.DeleteFile(tmpfile)

	indexfile := file + ".bai"

	markdup_args := []string{"markdup", "-l", fmt.Sprintf("%d", comp),
		"-t", fmt.Sprintf("%d", nthreads), "--tmpdir=" + tmpDir,
		tmpfile, file}

	cmd_out, err := common.LogCommand(sambamba_bin, markdup_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error in sambamba markdup: %s\n", string(cmd_out))
		return nil, err
	}

	return []string{indexfile}, nil
}
