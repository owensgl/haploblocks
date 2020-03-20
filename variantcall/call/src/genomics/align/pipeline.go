package main

import "fmt"
import "time"
import "log"
import "os"
import "path/filepath"
import "genomics/common"

type PipelineOp func(string, string, string, bool, int, string) ([]string, error)
type IndexOp func(string) (string, error)
type SortOp func(string, int, string) error

func lossyPipeline(s1, s2 string, output_dir string, finalize bool, compression int, tmpDir string) ([]string, error) {

	var files []string
	var runSort SortOp
	var runIndex IndexOp

	runSort = samSort
	runIndex = samIndex
	if sambam == true {
		runSort = sambamSort
		runIndex = sambamIndex
	}

	path_prefix := filepath.Join(output_dir, sample)

	// Trim low quality samples
	log.Println("Trimming samples ...")
	trim_start := time.Now()
	r1p, r2p, r1up, r2up, err := trimmomatic(s1, s2, output_dir, trim_prefix, false)
	if err != nil {
		log.Printf("Cannot trim samples: %s\n", err.Error())
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(trim_start))

	// Align paired and unpaired files from trimmomatic
	r1r2_paired := path_prefix + ".paired.sam"
	r1_unpaired := path_prefix + ".R1.unpaired.sam"
	r2_unpaired := path_prefix + ".R2.unpaired.sam"

	log.Println("Aligning samples ...")
	align_start := time.Now()

	aerr := ngmAlign(r1p, r2p, r1r2_paired, true, true)

	if aerr == nil {
		if common.NonZeroFileExists(r1up) {
			aerr = ngmAlign(r1up, "", r1_unpaired, false, true)
		} else {
			log.Printf("\tFile %s is empty, skipping ngmAlign...\n", r1up)
		}
	}

	if aerr == nil {
		if common.NonZeroFileExists(r2up) {
			aerr = ngmAlign(r2up, "", r2_unpaired, false, true)
		} else {
			log.Printf("\tFile %s is empty, skipping ngmAlign...\n", r2up)
		}
	}

	if aerr != nil {
		log.Printf("Cannot align samples: %s\n", aerr.Error())
		return nil, fmt.Errorf("NextGenMap error")
	}
	log.Printf("Time to align samples: %s\n", time.Since(align_start))

	// Convert to BAM
	r1r2_paired_bam := path_prefix + ".paired.bam"
	r1_unpaired_bam := path_prefix + ".R1.unpaired.bam"
	r2_unpaired_bam := path_prefix + ".R2.unpaired.bam"

	log.Printf("Converting to BAM files ...")
	bam_start := time.Now()

	err1 := samToBam(r1r2_paired, r1r2_paired_bam, true)
	err2 := samToBam(r1_unpaired, r1_unpaired_bam, true)
	err3 := samToBam(r2_unpaired, r2_unpaired_bam, true)

	if err1 != nil || err2 != nil || err3 != nil {
		log.Printf("Cannot covert SAM(%s): %s\n", r2_unpaired, err.Error())
		return nil, fmt.Errorf("samtools to bam error")
	}
	log.Printf("\t\t\t done: %s\n", time.Since(bam_start))

	// Concatenate BAMs
	log.Printf("Concatenating BAM files ...")
	cat_name := path_prefix + ".bam"
	cat_start := time.Now()
	cat_args := []string{"cat", "-o", cat_name, r1r2_paired_bam,
		r1_unpaired_bam, r2_unpaired_bam}

	if cmd_out, err := common.LogCommand("samtools", cat_args...).CombinedOutput(); err != nil {
		log.Printf("Cannot concat SAMs: %s\n", string(cmd_out))
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(cat_start))

	common.DeleteFile(r1r2_paired_bam)
	common.DeleteFile(r1_unpaired_bam)
	common.DeleteFile(r2_unpaired_bam)

	// If this is the final step, set appropriate compression level
	comp := 0
	if cull_bed == "" || finalize == true {
		comp = compression
	}

	// Sort aligned BAM
	log.Printf("Sorting BAM file ...")
	sort_start := time.Now()
	if err := runSort(cat_name, comp, tmpDir); err != nil {
		log.Printf("Cannot sort BAMs: %s\n", err.Error())
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(sort_start))

	// Remove non-repetitive regions
	if cull_bed != "" {

		// Index the BAM file
		log.Printf("Indexing BAM file ...")
		index_start := time.Now()
		index_file, err := runIndex(cat_name)
		if err != nil {
			log.Printf("Cannot index BAM: %s\n", err.Error())
			return nil, err
		}
		log.Printf("\t\t\t done: %s\n", time.Since(index_start))

		culled_name := path_prefix + ".culled.bam"
		cull_start := time.Now()
		if err = cullNonrepetitiveRegions(cat_name, cull_bed, culled_name); err != nil {
			log.Printf("Cannot cull BAM: %s\n", err.Error())
			return nil, err
		}
		log.Printf("\t\t\t done: %s\n", time.Since(cull_start))

		common.DeleteFile(index_file)
		common.DeleteFile(cat_name)
		os.Rename(culled_name, cat_name)

		if finalize == true {
			comp = compression
		}

		log.Printf("Sorting BAM file ...")
		sort_start = time.Now()
		if err := runSort(cat_name, comp, tmpDir); err != nil {
			log.Printf("Cannot sort BAMs: %s\n", err.Error())
		}
		log.Printf("\t\t\t done: %s\n", time.Since(sort_start))
	}

	files = []string{cat_name}

	if finalize == true {
		log.Printf("Indexing BAM file ...")
		index_start := time.Now()
		index_file, err := runIndex(cat_name)
		if err != nil {
			log.Printf("Cannot index BAM: %s\n", err.Error())
			return nil, err
		}
		log.Printf("\t\t\t done: %s\n", time.Since(index_start))
		files = append(files, index_file)
	}

	return files, nil
}

func losslessPipeline(s1, s2 string, output_dir string, finalize bool, compression int, tmpDir string) ([]string, error) {

	var files []string
	var metrics string
	var err error

	path_prefix := filepath.Join(output_dir, sample)

	fastq := path_prefix + ".fastq"
	ubam := path_prefix + ".unaligned.bam"
	sam := path_prefix + ".aligned.sam"
	bam := path_prefix + ".aligned.bam"

	// Convert to unaligned bam
	/// TODO: Check with Greg if this naming matters or if all
	///       this metadata is overwritten with ngm anyway.
	log.Println("Convert to uBAM ...")
	ubam_start := time.Now()
	if err = picardFastqToSam(s1, s2, ubam, 0, true, sample,
		sample, sample, true); err != nil {

		log.Printf("Error converting Fastq to SAM: %s\n", err.Error())
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(ubam_start))

	/*
	   // Sort unaligned bam
	   if sambam == true {
	       log.Println("Sorting unaligned sample ...")
	       sort_start := time.Now()
	       if err = sambamSort(ubam, 0); err != nil {
	           log.Printf("Error sorting BAM: %s\n", err.Error())
	           return nil, err
	       }
	       log.Printf("\t\t\t done: %s\n", time.Since(sort_start))
	   }
	*/

	// Marking low-quality reads
	log.Println("Marking low-quality reads ...")
	illumina_start := time.Now()
	if _, metrics, err = picardMarkIlluminaAdapters(ubam, output_dir, 0); err != nil {
		log.Printf("Error marking Illumina adapters: %s\n", err.Error())
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(illumina_start))

	// Convert to fastq and give low-quality reads low score
	log.Println("Converting to fastq ...")
	conv_start := time.Now()
	if err = picardSamToFastq(ubam, fastq, false); err != nil {
		log.Printf("Error marking Illumina adapters: %s\n", err.Error())
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(conv_start))

	// Align
	log.Println("Aligning samples ...")
	align_start := time.Now()
	if err = ngmAlign(fastq, "", sam, true, true); err != nil {
		log.Printf("Error aligning the sample: %s\n", err.Error())
		return nil, err
	}
	log.Printf("Time to align samples: %s\n", time.Since(align_start))

	// Sort aligned sample
	log.Println("Sorting aligned sample ...")
	sort_start := time.Now()

	if sambam != true {
		if err = picardSortEx(sam, bam, 0); err != nil {
			log.Printf("Error sorting SAM: %s\n", err.Error())
			return nil, err
		}
	} else {
		if err = samToBam(sam, bam, true); err != nil {
			log.Printf("Error converting SAM: %s\n", err.Error())
			return nil, err
		}

		if err = sambamSortPicard(bam, 0, tmpDir); err != nil {
			log.Printf("Error sorting BAM: %s\n", err.Error())
			return nil, err
		}
	}
	log.Printf("\t\t\t done: %s\n", time.Since(sort_start))

	// Merge aligned and unaligned samples
	if finalize == false {
		compression = 0
	}

	merged := path_prefix + ".bam"
	log.Println("Merging aligned and unaligned samples ...")
	merge_start := time.Now()
	files, err = picardMergeBamAlignment(bam, ubam, merged, ref_genome, compression,
		finalize, finalize, true)
	if err != nil {
		log.Printf("Error merging SAMs: %s\n", err.Error())
		return nil, err
	}
	log.Printf("\t\t\t done: %s\n", time.Since(merge_start))

	// Append the illumina metrics to final file list
	files = append(files, metrics)

	return files, nil
}
