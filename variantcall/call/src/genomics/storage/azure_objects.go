package storage

import "io"
import "os"

import "github.com/Azure/azure-sdk-for-go/storage"

/* Azure blobs are the same as S3 and GCP objects. */
func (as *AzureStorage) ListObjects(bname string) ([]string, error) {

    var objects []string

    bucket := as.client.GetContainerReference(bname)

    objs, err := bucket.ListBlobs(storage.ListBlobsParameters{})
    if err != nil {
        return nil, err
    }

    for _, b := range(objs.Blobs) {
        objects = append(objects, b.Name)
    }

    return objects, nil
}

func (as *AzureStorage) CopyTo(src, bucket, dst string) (error) {

    // Upload in 2MB chunks
    const bufsize = 1024 * 1024 * 2

    fd, err := os.Open(src)
    if err != nil {
        return err
    }
    defer fd.Close()

    blob := as.client.GetContainerReference(bucket).GetBlobReference(dst)
    if err := blob.PutAppendBlob(nil); err != nil {
        return err
    }

    buf := make([]byte, bufsize)
    eof := false

    for eof != true {
        n, err := fd.Read(buf)
        if err != nil {
            return err
        }

        if n < bufsize {
            buf = buf[:n]
            eof = true
        }

        if err := blob.AppendBlock(buf, nil); err != nil {
            return err
        }
    }

    return nil
}

func (as *AzureStorage) CopyFrom(bucket, src, dst string) (error) {

    blob := as.client.GetContainerReference(bucket).GetBlobReference(src)

    reader, err := blob.Get(nil)
    if err != nil {
        return err
    }
    defer reader.Close()

    fd, err := os.Create(dst)
    if err != nil {
        return err
    }

    if _, err = io.Copy(fd, reader); err != nil {
        return err
    }
    fd.Sync()

    if err = fd.Close(); err != nil {
        return err
    }

    return nil
}

func (as *AzureStorage) DeleteObject(bucket, obj string) (error) {

    blob := as.client.GetContainerReference(bucket).GetBlobReference(obj)
    _, err := blob.DeleteIfExists(&storage.DeleteBlobOptions{})
    return err
}
