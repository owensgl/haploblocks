package storage

import "github.com/Azure/azure-sdk-for-go/storage"

/* Azure storage containers are the same as S3 and GCP
 * buckets. Stick to the buckets terminology to avoid
 * confusion with compute containers.
 */
func (as *AzureStorage) ListBuckets() ([]string, error) {

    var buckets []string

    resp, err := as.client.ListContainers(storage.ListContainersParameters{})
    if err != nil {
        return nil, err
    }

    for _, container := range(resp.Containers) {
        buckets = append(buckets, container.Name)
    }

    return buckets, nil
}

func (as *AzureStorage) BucketExists(bname string) (bool) {

    bucket := as.client.GetContainerReference(bname)
    exists, _ := bucket.Exists()

    return exists
}

func (as *AzureStorage) CreateBucket(bname string) (error) {

    opts := storage.CreateContainerOptions{Access: storage.ContainerAccessTypePrivate}

    bucket := as.client.GetContainerReference(bname)
    _, err := bucket.CreateIfNotExists(&opts)
    return err
}

func (as *AzureStorage) DeleteBucket(bname string) (error) {

    bucket := as.client.GetContainerReference(bname)
    return bucket.Delete(&storage.DeleteContainerOptions{})
}
