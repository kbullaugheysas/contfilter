package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"strings"
	"sync"
)

type BamScanner struct {
	LineNumber int
	filename   string
	stdin      bool
	scanner    *bufio.Scanner
	wg         sync.WaitGroup
	prev       string
	record     []string
	Closed     bool
}

func (s *BamScanner) OpenBam(bamfile string) error {
	s.filename = bamfile
	cmd := exec.Command("samtools", "view", bamfile)
	input, err := cmd.StdoutPipe()
	if err != nil {
		return fmt.Errorf("failed creating pipe: %v", err)
	}
	if err := cmd.Start(); err != nil {
		return fmt.Errorf("command failed to start: %v", err)
	}
	s.scanner = bufio.NewScanner(input)
	s.wg.Add(1)
	go func() {
		s.wg.Wait()

		if !s.stdin {
			if err := cmd.Wait(); err != nil {
				log.Fatal("wait failed: ", err)
			}
		}
	}()
	return nil
}

// Fast forward to the next record with read name `read`
func (s *BamScanner) Find(read string) ([]string, error) {
	for {
		// The end of the file may have been reached previously.
		if s.Closed {
			return nil, nil
		}
		record, err := s.Record()
		if err != nil {
			return nil, err
		}
		// Or maybe the file is only now realized to be at the end.
		if s.Closed {
			return nil, nil
		}
		if record[0] == read {
			s.Ratchet()
			return record, nil
		}
		if strnum_cmp(record[0], read) < 0 {
			// Not far enough yet
			s.Ratchet()
		} else {
			// We didn't find the read before we reached one that is past what
			// we're looking for. We'll leave this one in the cache in case we
			// search for it next.
			return nil, nil
		}
	}
}

func (s *BamScanner) Record() ([]string, error) {
	if s.record != nil {
		return s.record, nil
	}
	s.Closed = !s.scanner.Scan()
	if err := s.scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner of %s errored: %v", s.filename, err)
	}
	if s.Closed {
		return nil, nil
	}
	line := strings.TrimSpace(s.scanner.Text())
	s.LineNumber++
	if len(line) == 0 {
		return nil, fmt.Errorf("empty BAM record")
	}
	s.record = strings.Split(line, "\t")
	if len(s.record) == 0 {
		return nil, fmt.Errorf("empty record at line %s", s.LineNumber)
	}
	read := s.record[0]
	if s.prev != "" {
		if strnum_cmp(s.prev, read) > 0 {
			return nil, fmt.Errorf("sorting order violated at line %d", s.LineNumber)
		}
	}
	s.prev = read
	return s.record, nil
}

func (s *BamScanner) Ratchet() {
	s.record = nil
}

func (s *BamScanner) Done() {
	s.wg.Done()
}

func (s *BamScanner) OpenStdin() {
	s.filename = "stdin"
	s.stdin = true
	s.wg.Add(1)
	s.scanner = bufio.NewScanner(os.Stdin)
}

func ReadBamHeader(bamfile string) (string, error) {
	output, err := exec.Command("samtools", "view", "-H", bamfile).Output()
	if err != nil {
		return "", fmt.Errorf("failed to read header: %v", err)
	}
	return string(output), nil
}

type BamWriter struct {
	filename string
	wg       sync.WaitGroup
	fp       *os.File
}

func (w *BamWriter) Open(bamfile string) (io.WriteCloser, error) {
	w.filename = bamfile
	cmd := exec.Command("samtools", "view", "-b", "-o", bamfile, "-")
	fp, err := cmd.StdinPipe()
	if err != nil {
		return nil, fmt.Errorf("failed creating pipe: %v", err)
	}
	w.wg.Add(1)
	go func() {
		samOut, err := cmd.CombinedOutput()
		if len(samOut) > 0 {
			log.Println("samtools output:")
			log.Print(string(samOut))
		}
		if err != nil {
			log.Fatal("executing samtools for writing bam file failed: ", err)
		}
		w.wg.Done()
	}()
	return fp, nil
}

func (w *BamWriter) Wait() {
	w.wg.Wait()
}
