package main

import "fmt"
import "flag"
import "os"
import "io/ioutil"

import "genomics/common"
import "genomics/storage"

func usage() {
	fmt.Fprintf(os.Stderr, `Usage: cas [--type TYPE] --get DIGEST | --put FILE [CASDIR]+

  (cas == content addressible storage)

  Retrieves or stores a file based on its content digest,
  in the content directories in CASDIR. On a PUT, only
  the first CASDIR is updated. On a GET, all CASDIRs are
  looked up in order.

  Modes of operation:

   -get DIGEST   Retrieve full path of file with digest DIGEST.
                 Printed to stdout. "" is printed if the file is not
                 found.

   -put FILE     Store file FILE into the first CASDIR

  OPT is one of:
    -type TYPE    Type of digest, only 'md5' is supported.
    -creds CREDS  Load credentials from this file.
`)
}

func loadCaches(cacheDirs []string) []common.Cas {
	var caches []common.Cas
	for _, casdir := range cacheDirs {
		cache := common.DirCas(casdir)
		if cache == nil {
			fmt.Fprintf(os.Stderr, "Could not use cache dir %s\n", casdir)
			return nil
		}
		caches = append(caches, cache)
	}
	return caches
}

func mainWithCode() int {

	var digest, putfile string
	var digestType string
	var credsFile string
	var creds common.AuthenticationCreds

	flag.StringVar(&digest, "get", "", "Digest to retrieve")
	flag.StringVar(&putfile, "put", "", "Filename to store")
	flag.StringVar(&digestType, "type", "md5", "Type of digest (md5)")
	flag.StringVar(&credsFile, "creds", "", "Credentials file (yaml)")

	flag.Usage = usage
	flag.Parse()

	casdirs := flag.Args()

	if len(casdirs) < 1 {
		flag.PrintDefaults()
		return 1
	}

	if digest != "" && putfile != "" {
		fmt.Fprintln(os.Stderr, "you must specify either -get or -put, not both.")
		flag.Usage()
		return 1
	}

	if digest == "" && putfile == "" {
		fmt.Fprintln(os.Stderr, "you must specify either -get or -put")
		flag.Usage()
		return 1
	}

	creds, err := common.ParseCreds("", "", "", credsFile)
	if err != nil {
		fmt.Fprintln(os.Stderr, "cannot parse credentials file %s", credsFile)
		return 1
	}

	if putfile != "" {
		caches := loadCaches(casdirs[0:1])
		if caches == nil {
			return 1
		}

		tmpDir, err := ioutil.TempDir(".", ".cas-")
		if err != nil {
			fmt.Fprintf(os.Stderr, "Can't create working directory: %s", err.Error())
			return 1
		}
		storage.LocalDir = tmpDir
		defer os.Remove(tmpDir)

		file, err := storage.Get(putfile, creds, true)
		if err != nil {
			fmt.Fprintln(os.Stderr, "could not access file")
			return 1
		}
		defer os.Remove(file)

		stored, err := caches[0].Put(file)
		if err != nil {
			fmt.Fprintln(os.Stderr, "Error:", err)
			return 1
		}
		// print cached input
		fmt.Println(stored)
	} else {
		dTyp := common.DigestTypeOf(digestType)
		if dTyp == common.DIGEST_UNKNOWN {
			fmt.Fprintln(os.Stderr, "invalid digest type %s", digestType)
			return 1
		}

		caches := loadCaches(casdirs)
		if caches == nil {
			return 1
		}

		found := ""
		for _, cache := range caches {
			found = cache.Get(dTyp, digest)
			if found != "" {
				fmt.Println(found)
				return 0
			}
		}
		fmt.Println("")
		return 0
	}
	return 0
}

func main() {
	os.Exit(mainWithCode())
}
