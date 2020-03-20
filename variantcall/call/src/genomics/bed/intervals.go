package bed

import "fmt"
import "strings"
import "bytes"
import "strconv"
import "os"
import "bufio"
import "io/ioutil"

//import "log"
import "math"
import "path/filepath"
import "genomics/common"
import "genomics/storage"

func getIntervalNaming(srcPath, prefix string, intervals uint64) (string, string) {

	ext := filepath.Ext(srcPath)
	m := filepath.Base(srcPath)
	//basename without ext
	n := m[:len(m)-len(ext)]

	pre := prefix
	if pre != "" {
		pre += "_"
	}
	name := fmt.Sprintf("%s%s_intervals_%d_", pre, n, intervals)

	return name, ext
}

func writeIntervals(name string, buffer *bytes.Buffer) error {

	f, err := os.Create(name)
	if err != nil {
		return fmt.Errorf("Error opening output file: %s\n", err.Error())
	}

	if _, err := buffer.WriteTo(f); err != nil {
		return fmt.Errorf("Error writing out to file: %s\n", err.Error())
	}

	if err := f.Close(); err != nil {
		return fmt.Errorf("Error closing file: %s\n", err.Error())
	}

	return nil
}

// CreateIntervals:
//
//  create sub-bed files with the given prefix, with `interval_size` bp each.
//  regions with less than `min_bed_size` are not included in the output.
//
//
//  The margin helps prevent splitting of contigs across segment boundaries.
//  Each interval segment file will be    interval_size +/- margin_size bp
//
func CreateIntervals(srcPath string, dstPrefix string, interval_size uint64, margin_size uint64, min_bed_size uint64,
	op_path string, creds common.AuthenticationCreds) (int, string, string, error) {

	var buffer bytes.Buffer
	var tot_size uint64
	var line string

	name, ext := getIntervalNaming(srcPath, dstPrefix, interval_size)

	fd, err := os.Open(srcPath)
	if err != nil {
		return 0, "", "", fmt.Errorf("Error opening file (%s): %s\n", err.Error())
	}
	defer fd.Close()

	// BED files are <Chromosome  Start  End>
	count := 0
	tot_size = 0
	srcLineNo := 0
	line = ""
	scanner := bufio.NewScanner(fd)
	for {

		var start, end uint64
		var err error

		if line == "" {
			if scanner.Scan() == false {
				t := fmt.Sprintf("%s%d%s", name, count, ext)

				// write outstanding data, if any
				if buffer.Len() > 0 {
					//log.Printf("EOF and buffer not empty -- flushing... (%d bp)", tot_size)
					if err = writeIntervals(t, &buffer); err != nil {
						return 0, "", "", err
					}
					count++

					if op_path != "" {
						if err = storage.Put(op_path, t, creds, false); err != nil {
							return 0, "", "", err
						}
					}
				}
				break
			}
			line = scanner.Text()
			srcLineNo++
		}

		parts := strings.Fields(line)

		if len(parts) < 3 {
			return 0, "", "", fmt.Errorf("Corrupted BED file: %d columns in %s\n", len(parts), line)
		}

		if start, err = strconv.ParseUint(parts[1], 10, 64); err != nil {
			return 0, "", "", fmt.Errorf("Error parsing start index: %s\n", err.Error())
		}

		if end, err = strconv.ParseUint(parts[2], 10, 64); err != nil {
			return 0, "", "", fmt.Errorf("Error parsing end index: %s\n", err.Error())
		}

		size := end - start

		// Skip tiny regions
		if size < min_bed_size {
			line = ""
			continue
		}

		// If the total size is less than the interval, some regions are added
		if tot_size < (interval_size - margin_size) {
			var towrite string

			if (tot_size + size) < (interval_size + margin_size) {
				// If region fits within margin, add to buffer
				tot_size += size
				towrite = line
				line = ""

			} else {
				// If not in either margin, split interval
				space := interval_size - tot_size
				nend := start + space

				strstart := fmt.Sprintf("%d", start)
				strend := fmt.Sprintf("%d", end)
				strnend := fmt.Sprintf("%d", nend)

				tot_size += space
				linedata := strings.TrimPrefix(line, parts[0])

				//log.Printf("L%d line before: %s  (size=%d, space=%d, new_end=%d)", srcLineNo, line, size, space, nend)
				towrite = parts[0] + strings.Replace(linedata, strend, strnend, 1)
				line = parts[0] + strings.Replace(linedata, strstart, strnend, 1)
				//log.Printf("L%d line after:  %s", srcLineNo, line)
			}

			buffer.WriteString(towrite)
			buffer.WriteByte('\n')
		}

		// Check if is time to flush
		if tot_size >= (interval_size - margin_size) {
			t := fmt.Sprintf("%s%d%s", name, count, ext)
			if err = writeIntervals(t, &buffer); err != nil {
				return 0, "", "", err
			}

			if op_path != "" {
				if err = storage.Put(op_path, t, creds, false); err != nil {
					return 0, "", "", err
				}
			}

			//log.Printf("flushing. count=%d tot_size=%d interval_size=%d", count+1, tot_size, interval_size)
			tot_size = 0
			count++
			buffer.Reset()
		}
	}

	return count, name, ext, nil
}

func countBasePairs(file string) (uint64, error) {

	var count uint64

	fd, err := os.Open(file)
	if err != nil {
		return 0, fmt.Errorf("Error opening file (%s): %s\n", err.Error())
	}
	defer fd.Close()

	// BED files are <Chromosome  Start  End>
	count = 0
	scanner := bufio.NewScanner(fd)
	for scanner.Scan() {

		var start, end uint64

		line := scanner.Text()
		parts := strings.Fields(line)

		if len(parts) < 3 {
			return 0, fmt.Errorf("Corrupted BED file: %d columns in %s\n", len(parts), line)
		}

		if start, err = strconv.ParseUint(parts[1], 10, 64); err != nil {
			return 0, fmt.Errorf("Error parsing start index: %s\n", err.Error())
		}

		if end, err = strconv.ParseUint(parts[2], 10, 64); err != nil {
			return 0, fmt.Errorf("Error parsing end index: %s\n", err.Error())
		}

		count += end - start
	}

	return count, nil
}

func SplitIntervals(srcPath, dstPrefix string, segments int, min_bed_size uint64) ([]string, error) {

	var intervals []string

	total_bps, err := countBasePairs(srcPath)
	if err != nil {
		return nil, err
	}

	interval_size := total_bps / uint64(segments)
	margin_size := uint64(float64(interval_size) / math.Max(100.0, float64(segments)))

	n, name, ext, err := CreateIntervals(srcPath, dstPrefix, interval_size, margin_size, min_bed_size, "", common.AuthenticationCreds{})
	if err != nil {
		return nil, err
	}

	if n > segments+1 {
		return nil, fmt.Errorf("Too many segments created (n=%d)", n)
	}

	// Merge overflow into last file
	if n > segments {
		overflow := fmt.Sprintf("%s%d%s", name, n-1, ext)
		last := fmt.Sprintf("%s%d%s", name, segments-1, ext)

		data, err := ioutil.ReadFile(overflow)
		if err != nil {
			return nil, err
		}

		f, err := os.OpenFile(last, os.O_APPEND|os.O_WRONLY, 0644)
		if err != nil {
			return nil, err
		}

		if _, err := f.Write(data); err != nil {
			return nil, err
		}

		if err := f.Close(); err != nil {
			return nil, err
		}

		common.DeleteFile(overflow)
	}

	for i := 0; i < segments; i++ {
		f := fmt.Sprintf("%s%d%s", name, i, ext)

		if common.FileExists(f) == false {
			return nil, fmt.Errorf("Interval file missing: %s", err.Error())
		}

		intervals = append(intervals, f)
	}

	return intervals, nil
}
