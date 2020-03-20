package storage

import "flag"
import "os"
import "testing"
import "time"
import "math/rand"

import "genomics/common"

var config struct {
    creds common.AuthenticationCreds
    protocol string
}

func TestMain(m *testing.M) {
    var user, pass, project, protocol string

    flag.StringVar(&user, "u", "", "Project User")
    flag.StringVar(&pass, "p", "", "Project Password/Access Key")
    flag.StringVar(&project, "n", "", "Project Namespace")
    flag.StringVar(&protocol, "t", "", "Cloud Protocol (AWS [s3], Azure [azure], GCP [gs])")

    flag.Parse()

    config.creds = common.AuthenticationCreds{Username: user, Password: pass, Project: project}
    config.protocol = protocol

    rand.Seed(time.Now().UnixNano())
    os.Exit(m.Run())
}
