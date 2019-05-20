package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"time"
)

var sample string
var margin float64
var minsamplen int
var mincontlen int
var maxdist int
var limit int
var penalty float64
var output string
var filterErcc bool
var logfn string

func init() {
	log.SetFlags(0)
	flag.StringVar(&sample, "sample", "", "BAM file of the sample you want to filter (sorted by name, required)")
	flag.Float64Var(&margin, "margin", 1.0, "how much better sample needs to be matched")
	flag.IntVar(&minsamplen, "min-samp-len", 60, "min length for a sample match")
	flag.IntVar(&mincontlen, "min-cont-len", 70, "min length for a contamination match")
	flag.IntVar(&maxdist, "max-edit-dist", 5, "max edit distance for a sample match")
	flag.IntVar(&limit, "limit", 0, "limit the number of sample reads considered (0 = no limit)")
	flag.Float64Var(&penalty, "edit-penalty", 2.0, "multiple for how to penalize edit distance")
	flag.StringVar(&output, "output", "", "output bam file (required)")
	flag.StringVar(&logfn, "log", "", "keep a record of what happens to each read")
	flag.BoolVar(&filterErcc, "ercc", false, "exclude ERCC mappings from sample before filtering")
	flag.Usage = func() {
		log.Println("usage: contfilter [options] cont1.bam cont2.bam")
		flag.PrintDefaults()
	}
}

func benchmark(start time.Time, label string) {
	elapsed := time.Since(start)
	log.Printf("%s took %s", label, elapsed)
}

func extract(row []string) (int, int, error) {
	if len(row) < 15 {
		return 0, 0, fmt.Errorf("too few fields")
	}
	match_len := len(row[9])
	edit_tag := row[14]
	if edit_tag[:5] != "nM:i:" {
		return 0, 0, fmt.Errorf("malformed edit distance tag: %s", edit_tag)
	}
	edit_dist, err := strconv.Atoi(edit_tag[5:])
	if err != nil {
		return 0, 0, fmt.Errorf("failed to parse edit dist: %s", edit_tag)
	}
	return match_len, edit_dist, nil
}

func main() {
	var kept_percent float64
	flag.Parse()
	contamination := flag.Args()
	startedAt := time.Now()

	if len(contamination) == 0 {
		log.Println("Must specify at least one contamination mapping BAM file")
		os.Exit(1)
	}

	if output == "" {
		log.Println("Must specify -output file")
		os.Exit(1)
	}

	scanner := BamScanner{}
	if sample == "" {
		scanner.OpenStdin()
	} else {
		if err := scanner.OpenBam(sample); err != nil {
			log.Fatal(err)
		}
	}

	reads_found := make([]int, len(contamination))
	reads_filtered := make([]int, len(contamination))
	contScanners := make([]BamScanner, len(contamination))
	rejected := make([]bool, len(contamination))
	found := make([]bool, len(contamination))

	for c := 0; c < len(contamination); c++ {
		if err := contScanners[c].OpenBam(contamination[c]); err != nil {
			log.Fatal(err)
		}
		reads_found[c] = 0
		reads_filtered[c] = 0
	}

	header, err := ReadBamHeader(sample)
	if err != nil {
		log.Fatal(err)
	}

	out := BamWriter{}
	outfp, err := out.Open(output)
	if err != nil {
		log.Fatal(err)
	}

	io.WriteString(outfp, header)

	var logfile *os.File
	if logfn != "" {
		log.Println("Opening log file", logfn)
		logfile, err = os.Create(logfn)
		if err != nil {
			log.Fatal(err)
		}
		defer logfile.Close()
		if logfile != nil {
			log.Println("Opened log file", logfn)
		}
	}

	reads_kept := 0
	total_reads := 0
	ercc := 0
	considered := 0
	too_short := 0
	too_diverged := 0

	err = func() error {
		defer scanner.Done()
		defer benchmark(startedAt, "processing")

		for {
			for c, _ := range contamination {
				rejected[c] = false
				found[c] = false
			}
			record, err := scanner.Record()
			if err != nil {
				return fmt.Errorf("failed to read from sample BAM: %v after %d lines", err, total_reads)
			}
			if scanner.Closed {
				return nil
			}
			scanner.Ratchet()

			if total_reads > 0 && total_reads%100000 == 0 {
				kept_percent = float64(reads_kept) / float64(considered) * 100
				log.Printf("considered %d out of %d so far, kept %0.1f%%\n", considered, total_reads, kept_percent)
			}
			if limit > 0 && limit == total_reads {
				return nil
			}
			total_reads++
			sample_len, sample_edit_dist, err := extract(record)
			if err != nil {
				return err
			}
			if sample_len < minsamplen {
				too_short++
				continue
			}
			if sample_edit_dist > maxdist {
				too_diverged++
				continue
			}
			sample_score := float64(sample_len) - float64(sample_edit_dist)*penalty
			read := record[0]
			if filterErcc && strings.Contains(record[2], "ERCC") {
				ercc++
				continue
			}
			considered++
			// Reads in the sample BAM will be rejected if either mate in any of the
			// contamination BAM files maps better than in the sampel BAM file.
			was_rejected := false
			for c := 0; c < len(contamination); c++ {
				for m := 0; m < 2; m++ {
					mate, err := contScanners[c].Find(read)
					if err != nil {
						log.Fatal(err)
					}
					if mate != nil {
						if !found[c] {
							found[c] = true
							reads_found[c]++
						}
						length, edit_dist, err := extract(mate)
						if err != nil {
							log.Fatalf("failed to read from %s: %v", contamination[c], err)
						}
						if length >= mincontlen {
							score := float64(length) - float64(edit_dist)*penalty
							if sample_score <= score+margin {
								if !rejected[c] {
									reads_filtered[c]++
									rejected[c] = true
									was_rejected = true
									if logfile != nil {
										entry := fmt.Sprintf("read %s with length %d and edit distance %d was rejected with score %0.1f because in %s it had a score of %0.1f with length %d and edit distance %d\n", read, sample_len, sample_edit_dist, sample_score, contamination[c], score, length, edit_dist)
										logfile.WriteString(entry)
									}
								}
							}
						}
					}
				}
			}
			if !was_rejected {
				// This read is okay, output it to the output BAM file.
				_, err := fmt.Fprintf(outfp, "%s\n", strings.Join(record, "\t"))
				if err != nil {
					return err
				}
				reads_kept++
				if logfile != nil {
					entry := fmt.Sprintf("Kept read %s with length %d and edit distance %d and score %0.1f\n",
						read, sample_len, sample_edit_dist, sample_score)
					logfile.WriteString(entry)
				}
			}
		}
	}()
	if err != nil {
		log.Fatal(err)
	}

	outfp.Close()
	out.Wait()

	if filterErcc {
		erccPerc := float64(ercc) / float64(total_reads) * 100
		log.Printf("filtered out %d ERCC reads (%0.1f%%) before comparing to contamination\n", ercc, erccPerc)
	}

	shortPerc := float64(too_short) / float64(total_reads) * 100
	log.Printf("filtered out %d reads (%0.1f%%) becase their alignment was too short\n", too_short, shortPerc)
	divergedPerc := float64(too_diverged) / float64(total_reads) * 100
	log.Printf("filtered out %d reads (%0.1f%%) becase they were too diverged\n", too_diverged, divergedPerc)

	for c, cont := range contamination {
		n := reads_filtered[c]
		perc := float64(n) / float64(considered) * 100
		found_perc := float64(reads_found[c]) / float64(considered) * 100
		log.Printf("found %d of %d reads in %s (%0.1f%%)\n", reads_found[c], considered, cont, found_perc)
		log.Printf("rejected %d of %d reads from %s (%0.1f%%)\n", reads_filtered[c], considered, cont, perc)
	}

	kept_percent = float64(reads_kept) / float64(considered) * 100
	total_percent := float64(reads_kept) / float64(total_reads) * 100
	log.Printf("kept %d of %d reads (%0.1f%%), which is %0.1f%% of the total\n",
		reads_kept, considered, kept_percent, total_percent)

	log.Println("machine parsable stats:")
	stats := []int{total_reads, ercc, considered, reads_kept}
	stats = append(stats, reads_found...)
	stats = append(stats, reads_filtered...)
	fmt.Print("stats")
	for _, s := range stats {
		fmt.Printf("\t%d", s)
	}
	fmt.Print("\n")
}
