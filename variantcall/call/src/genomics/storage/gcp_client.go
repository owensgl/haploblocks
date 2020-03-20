package storage

//import "io/ioutil"
//import "encoding/json"

import "golang.org/x/net/context"
import "cloud.google.com/go/storage"
//import "google.golang.org/api/option"

type GCPStorage struct {
    client *storage.Client
    project string
}

func (gs *GCPStorage) Create(project string) (error) {

    //var token interface{}

    client, err := storage.NewClient(context.Background())
                        //,
                        //option.WithServiceAccountFile(creds))
    if err != nil {
        return err
    }

/*
    data, err := ioutil.ReadFile(creds)
    if err != nil {
        return err
    }

    if err = json.Unmarshal(data, &token); err != nil {
        return err
    }
    project := token.(map[string]interface{})["project_id"].(string)
*/

    gs.client = client
    gs.project = project

    return nil
}
