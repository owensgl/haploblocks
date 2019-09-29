package storage

import "fmt"
import "strings"
import "log"
import "os"
import "path/filepath"

import "genomics/common"

func putLocalFile(path, file string, copy_file bool) error {

	dir := strings.TrimRight(strings.TrimLeft(path, "file:"), "/")
	dst := filepath.Join(dir, filepath.Base(file))

	if copy_file == false {

		// Move file
		if err := os.Rename(file, dst); err != nil {
			log.Printf("Error moving file: %s\n", err.Error())
			return err
		}

	} else {
		return common.CopyFile(file, dst)
	}

	return nil
}

func putCloudFile(path, file string, creds common.AuthenticationCreds) error {

	var store CloudStorage
	var err error

	proto, bucket, dir, _ := ParseCloudPath(path)

	if store, err = CreateCloudStore(proto, creds); err != nil {
		return err
	}

	dst := ""
	if dir != "" {
		dst = strings.TrimRight(dir, "/") + "/"
	}
	dst += filepath.Base(file)

	if err := store.CopyTo(file, bucket, dst); err != nil {
		log.Printf("Error uploading file (%s to %s): %s\n", file, path, err.Error())
		return err
	}

	return nil
}

func putFile(path, file string, creds common.AuthenticationCreds, copy_file bool) error {

	var err error

	backend := ParseBackend(strings.Split(path, ":")[0])
	if backend == LOCAL {
		err = putLocalFile(path, file, copy_file)
	} else if IsCloudBackend(backend) == true {
		err = putCloudFile(path, file, creds)
	} else {
		err = fmt.Errorf("%d is unsupported backend", backend)
	}

	return err
}
