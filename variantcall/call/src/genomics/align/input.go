package main

import "fmt"
import "strings"
import "unicode"
import "log"
import "time"
import "bufio"
import "os"
import "path/filepath"
import "net/url"
import "regexp"
import "genomics/common"
import "genomics/storage"

// reads entry named entryName from md5sum file digestFile
// returns "" if not found
func getMD5Entry(digestFile string, entryName string) string {
	var words []string
	var text string
	digestFd, err := os.Open(digestFile)
	if err != nil {
		log.Printf("can't open digestFile: %s", err.Error())
		return ""
	}
	defer digestFd.Close()

	scanner := bufio.NewScanner(digestFd)
	delim := regexp.MustCompile(`\s+`)
	expectedDigest := ""

	for scanner.Scan() {
		text = strings.TrimLeftFunc(scanner.Text(), unicode.IsSpace)
		words = delim.Split(text, 2)
		if len(words) < 2 {
			continue
		}
		if entryName == strings.TrimRight(words[1], "\r\n") {
			expectedDigest = words[0]
			break
		}
	}
	return expectedDigest
}

func verifyFileUsingDigestFile(digestFile string, entryName string, dataFile string) error {

	expectedDigest := getMD5Entry(digestFile, entryName)
	if expectedDigest == "" {
		return fmt.Errorf("Missing entry for name %s\n", entryName)
	}

	//log.Printf("Expecting file %s to have digest: MD5 %s", dataFile, expectedDigest)

	md5Sum, err := common.CalculateMD5Sum(dataFile)
	if err != nil {
		return err
	}

	// <digest>' '<fname>
	if md5Sum != expectedDigest {
		return fmt.Errorf("Checksum mismatch with Digest. expected %s but got %s\n", expectedDigest, md5Sum)
	}
	return nil
}

func verifyFileUsingHashFragment(dataFile string, hashFragment string) error {
	md5Sum, err := common.CalculateMD5Sum(dataFile)
	if err != nil {
		return err
	}

	if md5Sum != hashFragment {
		return fmt.Errorf("Checksum mismatch with hash fragment. expected %s but got %s\n", hashFragment, md5Sum)
	}
	return nil
}

func unzipFile(file string) (string, error) {

	cmd_out, err := common.LogCommand("gunzip", file).CombinedOutput()
	if err != nil {
		log.Printf("Error unzipping file: %s\n", string(cmd_out))
		return "", err
	}
	file = strings.TrimRight(file, ".gz")

	if common.FileExists(file) == false {
		log.Printf("Missing file(%s): %s\n", file, err.Error())
	}

	return file, nil
}

func getSample(sampleUrl, hashUrl string, creds common.AuthenticationCreds, unzip bool) (
	string, error) {

	var file string
	var hfile string
	var err error
	var expectedHash string
	var basename string
	var hashFragment string
	var sampleUrlWithoutFragmentHash string

	ts_start := time.Now()
	parsedUrl, err := url.Parse(sampleUrl)

	if err != nil {
		return "", err
	}

	basename = common.GetLastSubstring(parsedUrl.RequestURI(), "/")

	sampleUrlWithoutFragmentHash = parsedUrl.Scheme + "://" + parsedUrl.Host + parsedUrl.RequestURI()
	// TODO: maybe create a struct for the hash fragment
	_, hashFragment, err = common.ExtractHashFromUrlFragment(parsedUrl.Fragment)
	if err != nil {
		return "", err
	}
	if hashFragment != "" {
		// If hash is given in the fragment, use the hash from the fragment
		expectedHash = hashFragment
	} else if parsedUrl.Scheme == "hash" {
		expectedHash = parsedUrl.Host
	} else if hashUrl != "" {
		if hfile, err = storage.Get(hashUrl, creds, true); err != nil {
			log.Printf("Failed to get file (%s): %s\n", hashUrl, err.Error())
			return "", err
		}
		defer common.DeleteFile(hfile)

		// get the md5 hash from the url
		expectedHash = getMD5Entry(hfile, basename)

		if expectedHash == "" {
			log.Printf("Failed to determine MD5 digest from the hash url, %s, for %s", hashUrl, basename)
		}
	}

	if file, err = storage.Cached(common.DIGEST_MD5, expectedHash, creds, true); err != nil {
		log.Printf("Problem looking up content cache using the hash, %s , : %s", expectedHash, err.Error())
		file = ""
	}
	if file != "" {
		log.Printf("Retrieved content from cache. placed in %s", file)
	} else if parsedUrl.Scheme == "hash" {
		log.Printf("The cache folder doesn't have the hash %s", parsedUrl.EscapedPath())
	} else {
		if file, err = storage.Get(sampleUrlWithoutFragmentHash, creds, true); err != nil {
			log.Printf("Failed to get file (%s): %s\n", sample, err.Error())
			return "", err
		}

		if hashFragment != "" {
			log.Println("Verifying files using the hash fragment ...")
			if err = verifyFileUsingHashFragment(file, hashFragment); err != nil {
				log.Printf("Error verifying checksum with the hash fragment: %s\n", err.Error())
				return "", err
			}
		} else if hashUrl != "" {
			log.Println("Verifying files using the hash url ...")
			if err = verifyFileUsingDigestFile(hfile, basename, file); err != nil {
				log.Printf("Error verifying checksum with the hash url: %s\n", err.Error())
				return "", err
			}
		}
	}

	if unzip == true && strings.HasSuffix(basename, ".gz") {
		log.Println("Unzipping files ...")
		if file, err = unzipFile(file); err != nil {
			log.Printf("gzip error: %s\n", err.Error())
			return "", err
		}
	}
	log.Printf("Time to get sample: %s\n", time.Since(ts_start))

	return file, nil
}

func getReference(path string, creds common.AuthenticationCreds) (string, error) {

	ref_file, err := storage.Get(path, creds, false)
	if err != nil {
		log.Printf("Error fetching reference genome: %s\n", err.Error())
		return "", err
	}

	proto, _, _, fname := storage.ParseCloudPath(path)
	dir := filepath.Dir(path)
	ext := "." + common.GetLastSubstring(fname, ".")
	prefix := strings.TrimSuffix(fname, ext)

	if storage.IsCloudBackend(proto) == true {
		// If other files with the same prefix exist, get them too
		if files, err := storage.List(dir, creds); err == nil {

			for _, f := range files {
				if strings.HasPrefix(f, prefix) {
					exists := common.FileExists(f)
					if exists == false {
						log.Printf("Reference file missing, retrieving %s%s\n", dir, f)
						storage.Get(dir+f, creds, true)
					} else {
						log.Printf("Reference file exists, skipping: %s\n", f)
					}
				}
			}
		} else {
			log.Printf("Error listing directory contents: %s\n", err.Error())
		}
	}

	return ref_file, nil
}
