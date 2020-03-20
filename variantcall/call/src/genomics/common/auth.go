package common

import "log"
import "io/ioutil"
import "gopkg.in/yaml.v2"

type AuthenticationCreds struct {
	Username string
	Password string
	Project  string
}

//
// load an AuthenticationCreds object first from the file credsFile (yaml),
// then override settings using the other parameters (if any).
func ParseCreds(user string, pass string, project string, credsFile string) (AuthenticationCreds, error) {
	var settings AuthenticationCreds

	if credsFile != "" {

		log.Printf("Credentials read from file: %s", credsFile)
		credsData, err := ioutil.ReadFile(credsFile)

		if err != nil {
			return settings, err
		}

		//log.Printf("Creds Data: %v", credsData)

		err = yaml.Unmarshal(credsData, &settings)
		if err != nil {
			return settings, err
		}
	}

	if user != "" {
		settings.Username = user
	}
	if pass != "" {
		settings.Password = pass
	}
	if project != "" {
		settings.Project = project
	}
	return settings, nil
}
