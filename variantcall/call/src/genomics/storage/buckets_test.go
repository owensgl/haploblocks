package storage

import "fmt"
import "strings"
import "testing"

import "cloud.google.com/go/storage"

import "genomics/common"

func printBuckets(buckets []string, attrs []storage.BucketAttrs) {

    fmt.Println("Bucket List")

    for idx, bucket := range buckets {

        fmt.Printf("\t%s ", bucket)
        if attrs != nil {
            fmt.Printf(": %s and %s", attrs[idx].Location, attrs[idx].StorageClass)
        }
        fmt.Printf("\n")
    }
}

func checkBucketAttrs(buckets []string, attrs []storage.BucketAttrs,
            new_bucket string, zone string, storage_class string) bool {

    for idx, bucket := range buckets {

        if bucket == new_bucket &&
           zone == attrs[idx].Location &&
           storage_class == attrs[idx].StorageClass {

            return true
        }
    }

    return false
}

func TestBuckets(t *testing.T) {

    var store CloudStorage

    var test_bucket string
    var res bool

    test_bucket = common.GetRandomString(8)

    fmt.Printf("------ TESTING %s BUCKETS ------\n", strings.ToUpper(config.protocol))

    store, err := CreateCloudStore(ParseBackend(config.protocol), config.creds)
    if err != nil {
        t.Fatalf("Error creating storage client: %s\n", err.Error())
    }

    res = store.BucketExists(test_bucket)
    fmt.Printf("%s exists: %t\n", test_bucket, res)
    if res == true {
        t.FailNow()
    }

    fmt.Printf("Creating Bucket %s\n", test_bucket)
    if err := store.CreateBucket(test_bucket); err != nil {
        t.Fatalf("Error creating bucket %s: %s\n", test_bucket, err.Error())
    }

    res = store.BucketExists(test_bucket)
    if res == false {
        t.Fatalf("Cannot find created bucket\n")
    }

    if bkts, err := store.ListBuckets(); err == nil {
        printBuckets(bkts, nil)
    } else {
        t.Fatalf("Cannot retrieve bucket list: %s\n", err.Error())
    }

    //if bkts, attrs, err := store.ListBucketsEx(); err == nil {
    //    printBuckets(bkts, attrs)
    //} else {
    //    t.Fatalf("Cannot retrieve bucket metadata: %s\n", err.Error())
    //}

    fmt.Printf("Deleting Bucket %s\n", test_bucket)
    if err := store.DeleteBucket(test_bucket); err != nil {
        t.Fatalf("Error deleting bucket %s: %s\n", test_bucket, err.Error())
    }

    res = store.BucketExists(test_bucket)
    if res == true {
        t.Fatalf("Deleted bucket still exists\n")
    }

    if bkts, err := store.ListBuckets(); err == nil {
        printBuckets(bkts, nil)
    } else {
        t.Fatalf("Cannot retrieve bucket list: %s\n", err.Error())
    }

    // Skip second half since attributed are not cloud-vendor agnostic
    //test_bucket = common.GetRandomString(8)

    //zone := "US-WEST1"
    //storage_class := "REGIONAL"
    //fmt.Printf("Creating Bucket with extended attributes %s\n", test_bucket)
    //if err := store.CreateBucketEx(test_bucket,
    //                                zone, storage_class); err != nil {
    //    t.Fatalf("Error creating bucket %s: %s\n", test_bucket, err.Error())
    //}

    //if bkts, attrs, err := store.ListBucketsEx(); err == nil {
    //    if checkBucketAttrs(bkts, attrs, test_bucket,
    //                        zone, storage_class) == false {
    //        t.Fatal("Bucket does not match expected metadata\n")
    //    }
    //} else {
    //    t.Fatalf("Cannot retrieve bucket metadata: %s\n", err.Error())
    //}


    //fmt.Printf("Deleting Bucket %s\n", test_bucket)
    //if err := store.DeleteBucket(test_bucket); err != nil {
    //    t.Fatalf("Error deleting bucket %s: %s\n", test_bucket, err.Error())
    //}

    //res = store.BucketExists(test_bucket)
    //if res == true {
    //    t.Fatalf("Deleted bucket still exists\n")
    //}

    fmt.Println("------ DONE ------")
}
