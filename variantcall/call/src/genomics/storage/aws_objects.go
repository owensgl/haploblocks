package storage

import "os"

import "github.com/aws/aws-sdk-go/aws"
import "github.com/aws/aws-sdk-go/service/s3"
import "github.com/aws/aws-sdk-go/service/s3/s3manager"

func (s3c *S3Storage) ListObjects(bname string) ([]string, error) {

    var objects []string

    opts := &s3.ListObjectsInput{Bucket: aws.String(bname)}
    res, err := s3c.client.ListObjects(opts)
    if err != nil {
        return nil, err
    }

    for _, obj := range res.Contents {
        objects = append(objects, aws.StringValue(obj.Key))
    }

    return objects, nil
}

func (s3c *S3Storage) CopyTo(src, bucket, dst string) (error) {

    fd, err := os.Open(src)
    if err != nil {
        return err
    }
    defer fd.Close()

    uploader := s3manager.NewUploaderWithClient(s3c.client)
    upParams := &s3manager.UploadInput{
        Bucket: aws.String(bucket),
        Key:    aws.String(dst),
        Body:   fd,
    }

    _, err = uploader.Upload(upParams)
    return err
}

func (s3c *S3Storage) CopyFrom(bucket, src, dst string) (error) {

    fd, err := os.Create(dst)
    if err != nil {
        return err
    }

    downloader := s3manager.NewDownloaderWithClient(s3c.client)
    downParams := &s3.GetObjectInput{
        Bucket: aws.String(bucket),
        Key:    aws.String(src),
    }

    if _, err := downloader.Download(fd, downParams); err != nil {
        return err
    }

    if err = fd.Close(); err != nil {
        return err
    }

    return nil
}

func (s3c *S3Storage) DeleteObject(bucket, obj string) (error) {

    delobj := &s3.DeleteObjectInput{
        Bucket: aws.String(bucket),
        Key:    aws.String(obj),
    }
    _, err := s3c.client.DeleteObject(delobj)
    return err
}
