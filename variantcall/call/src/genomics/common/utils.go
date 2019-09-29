package common

import (
	"math/rand"
	"os/exec"
	"strings"
	"time"
	"log"
	"os"
)

type LogCmd struct {
	*exec.Cmd
}

func LogCommand(name string, args ...string) *LogCmd {
	cmd := exec.Command(name, args...)
	return &LogCmd{cmd}
}

func (cmd *LogCmd) CombinedOutput() ([]byte, error) {
	log.Printf("CMD Args:%v", cmd.Cmd.Args)
	cmd_dir := cmd.Cmd.Dir
	if cmd_dir == "" {
		cwd, err := os.Getwd()
		if err == nil {
			cmd_dir = cwd
		}
	}
	log.Printf("CMD Dir :%s", cmd_dir)
	cmd_start := time.Now()
	cmd_out, err := cmd.Cmd.CombinedOutput()
	log.Printf("CMD Took:%s", time.Since(cmd_start))
	if err != nil {
		log.Printf("CMD Err :%s", err.Error())
	}
	return cmd_out, err
}

func GetRandomString(n int) string {
	const letters = "abcdefghijklmnopqrstuvwxyz"

	b := make([]byte, n)

	for i := range b {
		b[i] = letters[rand.Int63()%int64(len(letters))]
	}

	return string(b)
}

func GetLastSubstring(str, sep string) string {
	parts := strings.Split(str, sep)
	return parts[len(parts)-1]
}

func CalculateMD5Sum(inputContent string) (string, error) {
	md5_args := []string{"--", inputContent}
	op, err := LogCommand("md5sum", md5_args...).CombinedOutput()
	if err != nil {
		return "", err
	}
	res := strings.Fields(string(op))[0]

	return res, nil
}
