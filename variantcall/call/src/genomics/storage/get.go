package storage

import "fmt"
import "strings"
import "log"
import "os"
import "io/ioutil"
import "path/filepath"
import URL "net/url"
import "genomics/common"

// local directory where copied/downloaded files end up
var LocalDir = "."

func copyToLocal(srcPath string, dstBasename string) (string, error) {
	localPath := filepath.Join(LocalDir, dstBasename)
	srcAbs, err := filepath.Abs(srcPath)
	if err != nil {
		log.Printf("Can't obtain absolute path for %s", srcPath)
		return "", err
	}
	err = os.Symlink(srcAbs, localPath)
	if err != nil {
		log.Printf("Error creating symlink: %s\n", err.Error())
		return "", err
	}
	return localPath, nil
}

func getLocalFile(file string, copy_file bool) (string, error) {

	path := strings.TrimPrefix(file, "file:")

	if copy_file == true {
		localName, err := copyToLocal(path, filepath.Base(path))
		if err != nil {
			return "", err
		}
		path = localName
	}

	return path, nil
}

func getHTTPFile(url string, creds common.AuthenticationCreds) (string, error) {

	var err error
	//var tmpfile *os.File

	credentials := fmt.Sprintf("j_username=%s&j_password=%s",
		URL.QueryEscape(creds.Username), URL.QueryEscape(creds.Password))

	// wget -q --no-cookies --no-check-certificate -c <fname> --post-file <passwd>
	wget_args := []string{"-q", "--no-cookies", "--no-check-certificate"}

	if creds.Username != "" && creds.Password != "" {
		// Create file to POST to wget
		passfile, err := ioutil.TempFile(LocalDir, ".post-form-")
		if err != nil {
			return "", err
		}
		if err := passfile.Close(); err != nil {
			return "", err
		}
		if err = ioutil.WriteFile(passfile.Name(), []byte(credentials), 0644); err != nil {
			log.Printf("Error creating credentials file: %s\n", err.Error())
			return "", err
		}
		defer common.DeleteFile(passfile.Name())
		wget_args = append(wget_args, "--post-file", passfile.Name())
	}

	basename := common.GetLastSubstring(url, "/")
	tmpfile, err := ioutil.TempFile(LocalDir, ".get-http-")

	if err != nil {
		return "", err
	}
	if err := tmpfile.Close(); err != nil {
		return "", err
	}

	wget_args = append(wget_args,
		"-c", url,
		"-O", tmpfile.Name())

	log.Printf("Downloading file %s\n", url)
	cmd_out, err := common.LogCommand("wget", wget_args...).CombinedOutput()
	if err != nil {
		log.Printf("Error downloading file: %s\n", string(cmd_out))
		return "", err
	}

	localpath := filepath.Join(LocalDir, basename)
	err = os.Rename(tmpfile.Name(), localpath)
	if err != nil {
		return "", err
	}
	log.Printf("File saved to: %s", localpath)
	return localpath, nil
}

func getCloudFile(path string, creds common.AuthenticationCreds) (string, error) {

	var store CloudStorage
	var err error

	proto, bucket, rfile, lfile := ParseCloudPath(path)

	if store, err = CreateCloudStore(proto, creds); err != nil {
		return "", err
	}

	localName := filepath.Join(LocalDir, lfile)
	if err = store.CopyFrom(bucket, rfile, localName); err != nil {
		log.Printf("Error downloading file (%s:%s): %s\n", bucket, rfile, err.Error())
		return "", err
	}

	return lfile, nil
}

func getCached(digestType common.DigestType, digest string, creds common.AuthenticationCreds, copy_file bool) (string, error) {
	found := ""

	if digestType == common.DIGEST_UNKNOWN {
		return "", nil
	}

	for _, cache := range CasCaches {
		found = cache.Get(digestType, digest)
		if found != "" {
			break
		}
	}

	// not cached
	if found == "" {
		return "", nil
	}

	if copy_file == true {
		localName, err := copyToLocal(found, filepath.Base(found))
		if err != nil {
			return "", err
		}
		found = localName
	}

	return found, nil

}

func getFile(path string, creds common.AuthenticationCreds, copy_file bool) (string, error) {

	var file string
	var err error

	backend := ParseBackend(strings.Split(path, ":")[0])
	if backend == LOCAL {
		file, err = getLocalFile(path, copy_file)
	} else if backend == HTTP {
		file, err = getHTTPFile(path, creds)
	} else if IsCloudBackend(backend) == true {
		file, err = getCloudFile(path, creds)
	} else {
		file = ""
		err = fmt.Errorf("%d is unsupported backend", backend)
	}

	if err == nil {
		file, err = filepath.Abs(file)
	}

	return file, err
}
