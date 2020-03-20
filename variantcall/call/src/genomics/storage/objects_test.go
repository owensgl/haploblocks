package storage

import "fmt"
import "strings"
import "os"
import "io/ioutil"
import "testing"

import "genomics/common"

func createTemp(filename string) (string, error) {

    data := common.GetRandomString(128)
    err := ioutil.WriteFile(filename, []byte(data), 0644)
    return data, err
}

func compareFile(filename, data string) (bool) {

    fdata, err := ioutil.ReadFile(filename)
    if err != nil {
        return false
    }

    return string(fdata) == data
}

func printObjects(objs []string) {
    for _, obj := range objs {
        fmt.Printf("\t%s\n", obj)
    }
}

func TestObjects(t *testing.T) {

    var store CloudStorage

    test_bucket := common.GetRandomString(8)
    test_object := common.GetRandomString(8)

    fmt.Printf("------ TESTING %s OBJECTS ------\n", strings.ToUpper(config.protocol))

    store, err := CreateCloudStore(ParseBackend(config.protocol), config.creds)
    if err != nil {
        t.Fatalf("Error creating storage client: %s\n", err.Error())
    }

    fmt.Printf("Create local file: %s\n", test_object)
    test_data, err := createTemp(test_object)
    if err != nil {
        t.Fatalf("Error creating storage client: %s\n", err.Error())
    }
    defer os.Remove(test_object)

    fmt.Printf("Create bucket: %s\n", test_bucket)
    if err = store.CreateBucket(test_bucket); err != nil {
        t.Fatalf("Error creating bucket %s: %s\n", test_bucket, err.Error())
    }
    defer store.DeleteBucket(test_bucket)

    fmt.Printf("Uploading file to %s/%s (bucket/object)\n", test_bucket, test_object)
    if err = store.CopyTo(test_object, test_bucket, test_object); err != nil {
        t.Fatalf("Error uploading file: %s\n", err.Error())
    }

    if objs, err := store.ListObjects(test_bucket); err == nil {
        fmt.Printf("Objects in %s\n", test_bucket)
        printObjects(objs)
    } else {
        t.Fatalf("Cannot retrieve object list: %s\n", err.Error())
    }

    test_dst := common.GetRandomString(8)
    fmt.Printf("Downloading object to %s\n", test_dst)
    if err = store.CopyFrom(test_bucket, test_object, test_dst); err != nil {
        t.Fatalf("Error downloading file: %s\n", err.Error())
    }
    defer os.Remove(test_dst)

    if res := compareFile(test_dst, test_data); res != true {
        t.Fatalf("Data mismatch in file\n")
    }

    fmt.Println("Deleting remote object")
    if err := store.DeleteObject(test_bucket, test_object); err != nil {
        t.Fatalf("Error deleting remote object: %s\n", err.Error())
    }

    fmt.Println("------ DONE ------")
}
