package storage

import "github.com/aws/aws-sdk-go/aws"
import "github.com/aws/aws-sdk-go/aws/awserr"
import "github.com/aws/aws-sdk-go/service/s3"

func (s3c *S3Storage) ListBuckets() ([]string, error) {

    var buckets []string

    res, err := s3c.client.ListBuckets(nil)
    if err != nil {
        return nil, err
    }

    for _, b := range res.Buckets {
        buckets = append(buckets, aws.StringValue(b.Name))
    }

    return buckets, nil
}

func (s3c *S3Storage) BucketExists(bname string) (bool) {

    input := &s3.HeadBucketInput{Bucket: aws.String(bname)}

    _, err := s3c.client.HeadBucket(input)
    if err != nil {
        aerr := err.(awserr.RequestFailure)
        // 404 Not Found
        if aerr.StatusCode() == 404 {
            return false
        }
    }

    return true
}

func (s3c *S3Storage) CreateBucket(bname string) (error) {

    opts := &s3.CreateBucketInput{Bucket: aws.String(bname)}

    _, err := s3c.client.CreateBucket(opts)
    return err
}

func (s3c *S3Storage) CreateBucketEx(bname, zone string) (error) {

    opts := &s3.CreateBucketInput{
        Bucket: aws.String(bname),
        CreateBucketConfiguration:
                &s3.CreateBucketConfiguration {
                    LocationConstraint: aws.String(zone),
                },
    }

    _, err := s3c.client.CreateBucket(opts)
    return err
}

func (s3c *S3Storage) DeleteBucket(bname string) (error) {
    opts := &s3.DeleteBucketInput {
        Bucket: aws.String(bname),
    }

    _, err := s3c.client.DeleteBucket(opts)
    return err
}
