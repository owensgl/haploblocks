package storage

import "fmt"
import "log"
import "strings"
import "genomics/common"

type StorageBackend uint8

const (
	LOCAL StorageBackend = iota + 1
	HTTP
	AWS_S3
	AZURE_STORAGE
	GCP_STORAGE
	UNKNOWN_STORAGE
)

// caching layer to shortcircuit downloads when content hashes
// are known
var CasCaches = []common.Cas{}

type CloudStorage interface {
	ListBuckets() ([]string, error)
	BucketExists(bname string) bool
	CreateBucket(bname string) error
	DeleteBucket(bname string) error

	ListObjects(bname string) ([]string, error)
	CopyTo(src, bucket, dst string) error
	CopyFrom(bucket, src, dst string) error
	DeleteObject(bucket, obj string) error
}

func CreateCloudStore(proto StorageBackend, creds common.AuthenticationCreds) (
	CloudStorage, error) {

	var cstore CloudStorage
	var err error

	cstore = nil
	err = nil

	if proto == AWS_S3 {
		var s3_client S3Storage

		err = s3_client.Create(creds.Username, creds.Password)
		cstore = &s3_client

	} else if proto == AZURE_STORAGE {
		var azure_client AzureStorage

		err = azure_client.Create(creds.Username, creds.Password)
		cstore = &azure_client

	} else if proto == GCP_STORAGE {
		var gcp_client GCPStorage

		err = gcp_client.Create(creds.Project)
		cstore = &gcp_client
	} else {
		err = fmt.Errorf("Unknown storage type")
	}

	return cstore, err
}

func IsCloudBackend(s StorageBackend) bool {
	return (s == AWS_S3 || s == AZURE_STORAGE || s == GCP_STORAGE)
}

// ParseBackend returns a StorageBackend enum based on its string
// representation. Recognized possible inputs are "file", "http", "s3", "azure", "gs".
// UNKNOWN_STORAGE is returned otherwise.
func ParseBackend(s string) StorageBackend {

	if s == "file" {
		return LOCAL
	} else if s == "http" || s == "https" {
		return HTTP
	} else if s == "s3" {
		return AWS_S3
	} else if s == "azure" {
		return AZURE_STORAGE
	} else if s == "gs" {
		return GCP_STORAGE
	}

	return UNKNOWN_STORAGE
}

func ParseCloudPath(object string) (StorageBackend, string, string, string) {
	log.Printf("ParseCloudPath %s", object)
	sep := ":"
	idx := strings.Index(object, sep)

	var proto StorageBackend
	var after_proto string
	if idx == -1 {
		proto = LOCAL
		after_proto = object
	} else {
		proto = ParseBackend(object[:idx])
		after_proto = object[idx+len(sep):]
	}

	parts := strings.Split(after_proto, "/")

	// assumes after_proto starts with '/'
	bucket := parts[1]
	fname := parts[len(parts)-1]

	if proto == LOCAL {
		return proto, bucket, after_proto, fname
	} else {
		return proto, bucket, strings.Join(parts[2:], "/"), fname
	}
}

func List(dir string, creds common.AuthenticationCreds) ([]string, error) {

	var store CloudStorage
	var err error

	proto, bucket, _, _ := ParseCloudPath(dir)
	if store, err = CreateCloudStore(proto, creds); err != nil {
		return nil, err
	}

	files, err := store.ListObjects(bucket)
	if err != nil {
		log.Printf("Error retrieving file list (%s): %s\n", bucket, err.Error())
		return nil, err
	}

	return files, nil
}

// GetCached obtains a local file representation of path if the
// content is available in the content caches.  If the content has
// not been cached before, then ("", nil) is returned.
func Cached(digestType common.DigestType, digest string, creds common.AuthenticationCreds, copy_file bool) (string, error) {
	return getCached(digestType, digest, creds, copy_file)
}

// Get retrieves a file based on its fully specified path and returns the name
// where it can be accessed from. path syntax is:
//
//  BACKEND:BACKEND-SPECIFIC-PATH
//
//  user and password are passed for authentication purposes over http(s).
//
// cloud_key is only for cloud backends. copy_file is for local files, and
// determines if the returned filename is a copy of the source.
//
func Get(path string, creds common.AuthenticationCreds, copy_file bool) (string, error) {
	return getFile(path, creds, copy_file)
}

// Write file designated by filename (output of Get), over target backend-specific path. path syntax is:
//
//  BACKEND:BACKEND-SPECIFIC-PATH
//
// - only local and cloud backends are supported.
// - target will be overwritten.
// - for local files, copy_file determines if the file should be copied (true) or moved (false).
func Put(path, file string, creds common.AuthenticationCreds, copy_file bool) error {
	return putFile(path, file, creds, copy_file)
}
