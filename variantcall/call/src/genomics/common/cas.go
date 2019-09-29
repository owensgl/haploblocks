//
// Content Addressible Storage
//
//   Cache and retrieve files based on their digests.
//
//   Only a directory backing and MD5 is implemented.
//
//
package common

import "os"
import "crypto/md5"
import "log"
import "io"
import "io/ioutil"
import "path/filepath"
import "errors"
import "fmt"
import "strings"

type DigestType int

const (
	DIGEST_MD5 DigestType = iota
	DIGEST_UNKNOWN
)

func (typ DigestType) String() string {
	names := [...]string{
		"MD5",
		"UNKNOWN",
	}
	if typ < DIGEST_MD5 || typ > DIGEST_UNKNOWN {
		return "UNKNOWN"
	}
	return names[typ]
}

func DigestTypeOf(str string) DigestType {
	str = strings.ToLower(str)
	switch str {
	case "md5":
		return DIGEST_MD5
	default:
		return DIGEST_UNKNOWN
	}
}

func ExtractHashFromUrlFragment(fragment string) (DigestType, string, error) {
	if fragment == "" {
		return DIGEST_UNKNOWN, "", nil
	}
	if !strings.Contains(fragment, "=") {
		return DIGEST_UNKNOWN, "", errors.New(fragment + " is not in a proper format ")
	}

	digestTypeString := fragment[:strings.Index(fragment, "=")]
	digestType := DigestTypeOf(digestTypeString)

	if digestType == DIGEST_UNKNOWN {
		return DIGEST_UNKNOWN, "", errors.New("System doesn't recognize the hash type " + digestTypeString)
	}

	hashFragment := fragment[strings.Index(fragment, "=")+1:]

	return digestType, hashFragment, nil
}

type Cas interface {
	Put(source string) (string, error)
	Get(typ DigestType, digest string) string
}


type dirCas struct {
	cdir string
}

// Obtains a Cas from an existing directory
func DirCas(dir string) Cas {
	dirBase := dirCas{cdir: dir}
	return &dirBase
}

func (cache dirCas) String() string {
	return fmt.Sprintf("dir:%s", cache.cdir)
}

func (cache dirCas) _createDir() error {
	err := os.MkdirAll(cache.cdir, os.ModeDir|os.ModePerm)
	return err
}

//
// Stores the given file identified by fpath, and returns a path
// reference to a file that holds the content.
//
// The basename of the file returned may not be the same as fpath's.
//
// Returns the pathname, and an error in case it could not be stored.
//
func (cache dirCas) Put(fpath string) (string, error) {
	f, err := os.Open(fpath)
	defer f.Close()

	if err != nil {
		return "", err
	}

	dst_basename := filepath.Base(fpath)
	dst_fullpath := filepath.Join(cache.cdir, dst_basename)

	tmp_dst, err := ioutil.TempFile(cache.cdir, ".cas-write-")

	if err != nil {
		err = cache._createDir()
		if err != nil {
			tmp_dst, err = ioutil.TempFile(cache.cdir, ".cas-write-")
		}
	}

	if err != nil {
		return "", err
	}

	defer func() {
		os.Remove(tmp_dst.Name())
		os.Remove(tmp_dst.Name() + ".md5")
		os.Remove(tmp_dst.Name() + ".md5link")
	}()

	// perform the copy -- and digest at the same time
	h := md5.New()
	// no internal buffering
	tee := io.TeeReader(f, h)
	if _, err := io.Copy(tmp_dst, tee); err != nil {
		return "", err
	}
	tmp_dst.Close()

	md5_digest := fmt.Sprintf("%x", h.Sum(nil))
	log.Printf("%s %s", md5_digest, dst_basename)

	cached_file := cache.Get(DIGEST_MD5, md5_digest)
	if cached_file != "" {
		log.Printf("already cached.")
		return cached_file, nil
	}

	// FIXME lame.
	// allow two different filenames with same content
	if _, err := os.Stat(dst_fullpath); err == nil {
		return "", errors.New(fmt.Sprintf("destination file %s already exists.", dst_fullpath))
	}

	//
	// Add new entry
	//

	digest_link := filepath.Join(cache.cdir, "MD5_"+md5_digest)
	if err = os.Symlink(dst_basename, tmp_dst.Name()+".md5link"); err != nil {
		return "", err
	}

	// file_foo
	// MD5_XXXXXXX ---> file_foo
	// file_foo.md5

	digest_data := []byte(fmt.Sprintf("%s\t%s\n", md5_digest, dst_basename))
	if err = ioutil.WriteFile(tmp_dst.Name()+".md5", digest_data, 0644); err != nil {
		return "", err
	}

	// rename temp files to their final location
	if err = os.Rename(tmp_dst.Name()+".md5link", digest_link); err != nil {
		return "", err
	}
	if err = os.Rename(tmp_dst.Name()+".md5", dst_fullpath+".md5"); err != nil {
		return "", err
	}
	// do this one last to seal the deal
	if err = os.Rename(tmp_dst.Name(), dst_fullpath); err != nil {
		return "", err
	}
	return dst_fullpath, nil
}

//
//  Returns the full path to file with the given content digest, or ""
//  if none could be found. The basename of the file returned may not
//  be the one passed to Put() -- do not depend on filenames being
//  preserved.
//
func (cache dirCas) Get(typ DigestType, digest string) string {
	if typ != DIGEST_MD5 {
		return ""
	}
	digest_link := filepath.Join(cache.cdir, typ.String()+"_"+digest)
	link_dest, err := os.Readlink(digest_link)
	digest = strings.TrimSpace(digest)

	if err != nil {
		return ""
	}

	full_path := ""
	if filepath.IsAbs(link_dest) {
		full_path = link_dest
	} else {
		full_path = filepath.Join(cache.cdir, link_dest)
	}

	_, err = os.Stat(full_path)
	if err != nil {
		return ""
	}
	return full_path
}
