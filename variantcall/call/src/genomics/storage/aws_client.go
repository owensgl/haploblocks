package storage

import "github.com/aws/aws-sdk-go/aws"
import "github.com/aws/aws-sdk-go/aws/session"
import "github.com/aws/aws-sdk-go/service/s3"
import "github.com/aws/aws-sdk-go/aws/credentials"

type S3Storage struct {
    client *s3.S3
}

func (s3c *S3Storage) Create(user, pass string) (error) {

    sess, err := session.NewSession(&aws.Config{
        Region:      aws.String("us-west-2"),
        Credentials: credentials.NewStaticCredentials(user, pass, ""),
    })
    if err != nil {
        return err
    }

    s3c.client = s3.New(sess)

    return nil
}
