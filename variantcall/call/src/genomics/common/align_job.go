package common

import "encoding/json"
import "io/ioutil"

type AlignSample struct {
    Name string `json:"name"`
    // Store in indexed as [read][file] in the form:
    // <file1, hash1>, <file2, hash2>
    Files [2][2] string `json:"locations"`
}

func (as *AlignSample) Init() {
    as.Name = ""
    as.Files[0][0] = ""
    as.Files[0][1] = ""
    as.Files[1][0] = ""
    as.Files[1][1] = ""
}

func LoadAlignJob(file string) (map[string]AlignSample, error) {
    var data map[string]AlignSample

    marshalled, err := ioutil.ReadFile(file)
    if err != nil {
        return nil, err
    }

    if err = json.Unmarshal(marshalled, &data); err != nil {
        return nil, err
    }

    return data, err
}

func StoreAlignJob(data map[string]AlignSample, file string) (error) {

    if marshalled, err := json.Marshal(data); err == nil {
        if err = ioutil.WriteFile(file, []byte(marshalled), 0644); err != nil {
            return err
        }
    } else {
        return err
    }

    return nil
}

func StoreSingleAlignJob(entry AlignSample, file string) (error) {

    data := make(map[string]AlignSample)
    data[entry.Name] = entry

    return StoreAlignJob(data, file)
}
