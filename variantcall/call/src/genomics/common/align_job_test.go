package common

import "fmt"
import "os"
import "strconv"
import "testing"

func checkAlignEntry(entry AlignSample, idx int) (error) {

    mismatch := fmt.Errorf("Mismatch in stored value\n")

    if v, err := strconv.Atoi(entry.Name); err == nil {
        if v != idx {
            return mismatch
        }
    } else {
        return err
    }

    for i := 0; i < 2; i++ {
        for j := 0; j < 2; j++ {
            if v, err := strconv.Atoi(entry.Files[i][j]); err == nil {
                if v != (idx * 100) + (2 * i) + j {
                    return mismatch
                }
            } else {
                return err
            }
        }
    }

    return nil
}

func TestAlignSamples(t *testing.T) {

    file := "TEST"
    data := make(map[string]AlignSample)
    count := 5
    tidx := 4

    for i := 0; i < count; i++ {
        var entry AlignSample

        entry.Name = strconv.Itoa(i)
        entry.Files[0][0] = strconv.Itoa((i * 100))
        entry.Files[0][1] = strconv.Itoa((i * 100) + 1)
        entry.Files[1][0] = strconv.Itoa((i * 100) + 2)
        entry.Files[1][1] = strconv.Itoa((i * 100) + 3)

        data[entry.Name] = entry
    }

    StoreAlignJob(data, file)
    defer os.Remove(file)

    stored, err := LoadAlignJob(file)
    if err != nil {
        t.Fatalf("Cannot retrieve align job: %s\n", err.Error())
    }

    for i := 0; i < count; i++ {
        key := strconv.Itoa(i)
        entry, res := stored[key]

        if res == false {
            t.Fatalf("Cannot retrieve stored entry: %s\n", err.Error())
        }

        if err := checkAlignEntry(entry, i); err != nil {
            t.Fatalf("Error in retrieved job: %s\n", err.Error())
        }
    }

    os.Remove(file)
    entry, _ := data[strconv.Itoa(tidx)]
    if err = StoreSingleAlignJob(entry, file); err != nil {
        t.Fatalf("Error in saving single entry: %s\n", err.Error())
    }

    ndata, err := LoadAlignJob(file)
    if err != nil {
        t.Fatalf("Error reading back single entry: %s\n")
    }

    if len(ndata) != 1 {
        t.Fatalf("More than a single entry read back\n")
    }

    if nentry, res := ndata[strconv.Itoa(tidx)]; res == true {
        if err = checkAlignEntry(nentry, tidx); err != nil {
            t.Fatalf("Error in single entry data: %s\n", err.Error())
        }
    } else {
        t.Fatalf("%d entry missing in single entry read back\n", tidx)
    }
}
