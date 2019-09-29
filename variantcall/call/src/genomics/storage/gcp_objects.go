package storage

import "io"
import "os"

import "golang.org/x/net/context"
import "cloud.google.com/go/storage"
import "google.golang.org/api/iterator"

func (gs *GCPStorage) ListObjectsEx(bucket string) (
            []string, []storage.ObjectAttrs, error) {

    var objects []string
    var attrs []storage.ObjectAttrs

    it := gs.client.Bucket(bucket).Objects(context.Background(), nil)
    for {
        attr, err := it.Next()
        if err == iterator.Done {
            break
        }
        if err != nil {
            return nil, nil, err
        }
        objects = append(objects, attr.Name)
        attrs = append(attrs, *attr)
    }

    return objects, attrs, nil
}

func (gs *GCPStorage) ListObjects(bucket string) (
            []string, error) {

    objects, _, err := gs.ListObjectsEx(bucket)
    return objects, err
}

func (gs *GCPStorage) CopyTo(src, bucket, dst string) (error) {

    ctx := context.Background()
    fd, err := os.Open(src)
    if err != nil {
        return err
    }
    defer fd.Close()

    writer := gs.client.Bucket(bucket).Object(dst).NewWriter(ctx)

    if _, err = io.Copy(writer, fd); err != nil {
        return err
    }

    if err = writer.Close(); err != nil {
        return err
    }

    return nil
}

func (gs *GCPStorage) CopyFrom(bucket, src, dst string) (error) {

    ctx := context.Background()
    reader, err := gs.client.Bucket(bucket).Object(src).NewReader(ctx)
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

func (gs *GCPStorage) DeleteObject(bucket, obj string) (error) {

    return gs.client.Bucket(bucket).Object(obj).Delete(context.Background())
}
