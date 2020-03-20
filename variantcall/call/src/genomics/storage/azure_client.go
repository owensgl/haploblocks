package storage

import "github.com/Azure/azure-sdk-for-go/storage"

type AzureStorage struct {
    client storage.BlobStorageClient
}

func (as *AzureStorage) Create(user, pass string) (error) {

    client, err := storage.NewBasicClient(user, pass)
    if err != nil {
        return err
    }

    as.client = client.GetBlobService()

    return nil
}
