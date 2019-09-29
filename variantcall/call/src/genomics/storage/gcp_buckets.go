package storage

import "golang.org/x/net/context"
import "cloud.google.com/go/storage"
import "google.golang.org/api/iterator"

func (gs *GCPStorage) ListBucketsEx() (
            []string, []storage.BucketAttrs, error) {

    var buckets []string
    var attrs []storage.BucketAttrs

    it := gs.client.Buckets(context.Background(), gs.project)
    for {
        attr, err := it.Next()
        if err == iterator.Done {
            break
        }
        if err != nil {
            return nil, nil, err
        }
        buckets = append(buckets, attr.Name)
        attrs = append(attrs, *attr)
    }

    return buckets, attrs, nil
}

func (gs *GCPStorage) ListBuckets() ([]string, error) {

    buckets, _, err := gs.ListBucketsEx()
    return buckets, err
}

func (gs *GCPStorage) BucketExists(bname string) (bool) {

    bucket := gs.client.Bucket(bname)
    _, err := bucket.Attrs(context.Background())

    if err == storage.ErrBucketNotExist {
        return false
    }

    return true
}

func (gs *GCPStorage) CreateBucket(bname string) (error) {

    bucket := gs.client.Bucket(bname)
    return bucket.Create(context.Background(), gs.project, nil)
}

func (gs *GCPStorage) CreateBucketEx(bname, zone, class string) (error) {

    attrs := storage.BucketAttrs{Location: zone, StorageClass: class}

    bucket := gs.client.Bucket(bname)
    return bucket.Create(context.Background(), gs.project, &attrs)
}

func (gs *GCPStorage) DeleteBucket(bname string) (error) {

    bucket := gs.client.Bucket(bname)
    return bucket.Delete(context.Background())
}
