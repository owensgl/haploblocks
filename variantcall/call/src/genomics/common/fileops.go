package common

import "log"
import "io"
import "os"

func FileExists(filename string) bool {
	if _, err := os.Stat(filename); err != nil {
		return false
	}

	return true
}

// Returns true if disk file exists and has size > 0. false otherwise.
func NonZeroFileExists(filename string) bool {

	if info, err := os.Stat(filename); err == nil {
		if info.Size() > 0 {
			return true
		}
	}

	return false
}

func CopyFile(src, dst string) error {
	var fiSrc os.FileInfo
	fiDst, err := os.Stat(dst)
	if err == nil {
		fiSrc, err = os.Stat(src)
		if err != nil {
			// problem with src
			return err
		}
		if os.SameFile(fiSrc, fiDst) {
			// files are the same. trivial copy
			return nil
		}
		return err
	}

	s, err := os.Open(src)
	if err != nil {
		return err
	}

	defer s.Close()

	d, err := os.Create(dst)
	if err != nil {
		return err
	}

	if _, err = io.Copy(d, s); err != nil {
		return err
	}

	if err = d.Close(); err != nil {
		return err
	}

	return nil
}

func DeleteFile(filename string) {
	if err := os.Remove(filename); err != nil {
		log.Printf("Error deleting filename(%s): %s\n", filename, err.Error())
	}
}
