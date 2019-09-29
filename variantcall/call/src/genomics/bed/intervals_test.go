package bed

import "fmt"
import "bytes"
import "os"
import "io"
import "bufio"
import "io/ioutil"
import "testing"
import "strings"
import "strconv"
import "compress/gzip"
import "log"
import "genomics/common"
import "errors"

const test_file = "TEST.bed"
const bed_gz = "test.bed.gz"

func largeBed(t *testing.T) {
	output, err := os.Create(test_file)
	if err != nil {
		t.Fatalf("Cannot create file %s: %s\n", test_file, err.Error())
	}
	input, err := os.Open(bed_gz)
	if err != nil {
		t.Fatalf("Cannot open file %s: %s\n", err.Error())
	}

	gzread, err := gzip.NewReader(input)
	num_written, err := io.Copy(output, gzread)
	log.Printf("Wrote %d bytes to file %s", num_written, test_file)

	if err != nil {
		t.Fatalf("Cannot write file %s: %s\n", test_file, err.Error())
	}

	if err := input.Close(); err != nil {
		t.Fatalf("Error closing file: %s\n", bed_gz, err.Error())
	}
	if err := output.Close(); err != nil {
		t.Fatalf("Error closing file: %s\n", test_file, err.Error())
	}
}

func createBed(t *testing.T) {
	var buffer bytes.Buffer
	f, err := os.Create(test_file)
	if err != nil {
		t.Fatalf("Cannot retrieve bed file: %s\n", err.Error())
	}

	buffer.WriteString("Chr0\t1000\t1500\n")
	buffer.WriteString("Chr0\t2200\t2210\n")
	buffer.WriteString("Chr0\t2500\t2700\n")

	if _, err := buffer.WriteTo(f); err != nil {
		t.Fatalf("Error writing out to file: %s\n", err.Error())
	}

	if err := f.Close(); err != nil {
		t.Fatalf("Error closing file: %s\n", err.Error())
	}
}

func TestCreateIntervals(t *testing.T) {

	res := []string{"Chr0\t1000\t1300\n",
		"Chr0\t1300\t1500\nChr0\t2500\t2600\n",
		"Chr0\t2600\t2700\n"}

	createBed(t)

	count, name, ext, err := CreateIntervals(test_file, "", 300, 3, 20, "", common.AuthenticationCreds{})
	if err != nil {
		t.Fatalf("Error creating intervals: %s\n", err.Error())
	}

	if count != 3 {
		t.Fatalf("Incorrect interval parsing: expected %d, real %d\n", 3, count)
	}

	for i := 0; i < count; i++ {

		s := fmt.Sprintf("%s%d%s", name, i, ext)
		data, err := ioutil.ReadFile(s)
		if err != nil {
			t.Fatalf("Error reading file %s: %s\n", s, err.Error())
		}

		str := string(data)
		if str != res[i] {
			t.Fatalf("Data mismatch in %s: %s\n", s, str)
		}

		if err := os.Remove(s); err != nil {
			t.Fatalf("Error deleting files: %s\n", err.Error())
		}
	}

	if err := os.Remove(test_file); err != nil {
		t.Fatalf("Error deleting files: %s\n", err.Error())
	}
}

type ContigLine struct {
	chrom string
	start uint64
	end   uint64
}

func getNextContig(files []*bufio.Scanner) ([]*bufio.Scanner, ContigLine) {
	var start, end uint64
	var err error
	if len(files) == 0 {
		return nil, ContigLine{}
	}
	scanner := files[0]

	if scanner.Scan() == false {
		return getNextContig(files[1:])
	}
	line := scanner.Text()
	parts := strings.Fields(line)
	if len(parts) < 3 {
		return nil, ContigLine{}
	}
	if start, err = strconv.ParseUint(parts[1], 10, 64); err != nil {
		return nil, ContigLine{}
	}

	if end, err = strconv.ParseUint(parts[2], 10, 64); err != nil {
		return nil, ContigLine{}
	}

	return files, ContigLine{parts[0], start, end}
}

// checks if contig a starts with contig b. returns a - b
//
func startsWith(a ContigLine, b ContigLine) (ContigLine, error) {

	var mismatch = fmt.Errorf("expected %s %d %d, found %s %d %d",
		a.chrom, a.start, a.end, b.chrom, b.start, b.end)

	if a.chrom != b.chrom {
		return a, mismatch
	}

	if a.start != b.start {
		return a, mismatch
	}

	if a.end < b.end {
		return a, mismatch
	}

	if a.end == b.end {
		//log.Printf("MATCH, %s %d %d", a.chrom, a.start, a.end)
		return ContigLine{}, nil
	}

	return ContigLine{a.chrom, b.end, a.end}, nil
}

//
// goes through each contig in the bed file, in order.
// makes sure that the contig is represented in one of the segment.
// the function assumess that the segments are sorted.
//
func compareContigs(bedfile string, segments []string, t *testing.T) error {
	main_fd, err := os.Open(bedfile)
	if err != nil {
		return err
	}
	defer main_fd.Close()

	main_scan := bufio.NewScanner(main_fd)

	scan2 := []*bufio.Scanner{}
	split_fds := []*os.File{}
	defer func() {
		for _, fd := range split_fds {
			if fd != nil {
				fd.Close()
			}
		}
	}()

	for i, fil := range segments {
		fd, err := os.Open(fil)
		if err != nil {
			t.Errorf("can't open split file %d:%s", i, err.Error())
			return err
		}
		split_fds = append(split_fds, fd)
		scan2 = append(scan2, bufio.NewScanner(fd))
	}

	var contig1, contig2 ContigLine

	scan1 := []*bufio.Scanner{main_scan}
	for {
		scan1, contig1 = getNextContig(scan1)
		if scan1 == nil || contig1.chrom == "" {
			// done
			break
		}

		for {
			if scan2 == nil {
				return errors.New("could not match parent contig")
			}
			scan2, contig2 = getNextContig(scan2)
			if contig2.chrom == "" {
				return errors.New("EOF prematurely reached in split files")
			}
			contig1, err = startsWith(contig1, contig2)
			if err != nil {
				return err
			}
			if contig1.chrom == "" {
				// we've matched successfully
				break
			}
		}
	}
	return nil
}

func TestSplitIntervalsLarge(t *testing.T) {
	largeBed(t)
	defer os.Remove(test_file)

	total_bp, err := countBasePairs(test_file)
	if err != nil {
		t.Errorf("Could not count BP: %s", err.Error())
		return
	}

	var files []string
	test_splits := []int{1, 2, 3, 5, 7, 10, 15, 17, 32, 64, 113, 128, 256, 257, 1377}
	for _, num_splits := range test_splits {
		log.Printf("splitting large BED into %d intervals...", num_splits)

		// don't skip short intervals
		files, err = SplitIntervals(test_file, "", num_splits, 0)
		if err != nil {
			t.Errorf("problem splitting interval: %s", err.Error())
			goto cleanup
		}
		if len(files) != num_splits {
			t.Errorf("expected %d splits. got %d", num_splits, len(files))
			goto cleanup
		}

		// Check that we have the same count in sub files.
		split_bp := uint64(0)
		for i := 0; i < len(files); i++ {
			sub_bp, err := countBasePairs(files[i])
			if err != nil {
				t.Errorf("could not count bp in split file %d: %s", i, err.Error())
				goto cleanup
			}
			split_bp += sub_bp
		}

		err := compareContigs(test_file, files, t)
		if err != nil {
			t.Errorf("contigs mismatch.:%s", err.Error())
			goto cleanup
		}

		if split_bp != total_bp {
			t.Errorf("split files have %d bp, expected %d", split_bp, total_bp)
			goto cleanup
		}

		for i, fil := range files {
			os.Remove(fil)
			files[i] = ""
		}
	}
cleanup:
	// for _, fil := range files {
	// 	if fil != "" {
	// 		os.Remove(fil)
	// 	}
	// }
}

func TestSplitIntervals(t *testing.T) {

	res := []string{"Chr0\t1000\t1355\n",
		"Chr0\t1355\t1500\nChr0\t2500\t2700\n"}

	createBed(t)

	files, err := SplitIntervals(test_file, "", 2, 20)
	if err != nil {
		t.Fatalf("Error creating intervals: %s\n", err.Error())
	}

	if len(files) != 2 {
		t.Fatalf("Incorrect interval parsing: expected %d, real %d\n", 2, len(files))
	}

	for i := 0; i < len(files); i++ {

		data, err := ioutil.ReadFile(files[i])
		if err != nil {
			t.Fatalf("Error reading file %s: %s\n", files[i], err.Error())
		}

		str := string(data)
		if str != res[i] {
			t.Fatalf("Data mismatch in %s: %s\n", files[i], str)
		}

		if err := os.Remove(files[i]); err != nil {
			t.Fatalf("Error deleting files: %s\n", err.Error())
		}
	}

	if err := os.Remove(test_file); err != nil {
		t.Fatalf("Error deleting files: %s\n", err.Error())
	}
}
