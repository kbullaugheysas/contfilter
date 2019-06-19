package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"time"
)

type Args struct {
	Sample      string
	Margin      float64
	MinLength   int
	MaxDist     int
	Limit       int
	Penalty     float64
	Output      string
	Ercc        bool
	LogFilename string
	Verbose     bool
}

var args = Args{}
var logger *log.Logger

func init() {
	log.SetFlags(0)
	flag.StringVar(&args.Sample, "sample", "", "BAM file of the sample you want to filter (sorted by name, required)")
	flag.Float64Var(&args.Margin, "margin", 1.0, "how much better sample needs to be matched")
	flag.IntVar(&args.MinLength, "min-len", 60, "min length for an alignment")
	flag.IntVar(&args.MaxDist, "max-edit-dist", 5, "max edit distance for a sample match")
	flag.IntVar(&args.Limit, "limit", 0, "limit the number of sample reads considered (0 = no limit)")
	flag.Float64Var(&args.Penalty, "edit-penalty", 2.0, "multiple for how to penalize edit distance")
	flag.StringVar(&args.Output, "output", "", "output bam file (required)")
	flag.StringVar(&args.LogFilename, "log", "", "write parameters and stats to a log file")
	flag.BoolVar(&args.Verbose, "verbose", false, "keep a record of what happens to each read in the log (must give -log name)")
	flag.BoolVar(&args.Ercc, "ercc", false, "exclude ERCC mappings from sample before filtering")
	flag.Usage = func() {
		log.Println("usage: contfilter [options] cont1.bam cont2.bam")
		flag.PrintDefaults()
	}
}

func benchmark(start time.Time, label string) {
	elapsed := time.Since(start)
	logger.Printf("%s took %s", label, elapsed)
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

func OpenLogger() {
	if args.LogFilename == "" {
		logger = log.New(os.Stderr, "", 0)
	} else {
		logfile, err := os.Create(args.LogFilename)
		if err != nil {
			log.Fatal(err)
		}
		logger = log.New(logfile, "", 0)
	}
}

func LogArguments() {
	logger.Println("command:", strings.Join(os.Args, " "))
	blob, err := json.MarshalIndent(args, "", "    ")
	if err != nil {
		logger.Fatal("failed to marshal arguments")
	}
	logger.Println(string(blob))
}

func MatchesErcc(mate1, mate2 []string) bool {
	return args.Ercc &&
		(strings.Contains(mate1[2], "ERCC") || (mate2 != nil && strings.Contains(mate2[2], "ERCC")))
}

func main() {
	var kept_percent float64
	flag.Parse()
	contamination := flag.Args()
	startedAt := time.Now()

	if len(contamination) == 0 {
		logger.Println("must specify at least one contamination mapping BAM file")
		os.Exit(1)
	}

	if args.Output == "" {
		logger.Println("must specify -output file")
		os.Exit(1)
	}

	OpenLogger()
	LogArguments()

	scanner := BamScanner{}
	if args.Sample == "" {
		scanner.OpenStdin()
	} else {
		if err := scanner.OpenBam(args.Sample); err != nil {
			logger.Fatal(err)
		}
	}

	reads_found := make([]int, len(contamination))
	reads_filtered := make([]int, len(contamination))
	contScanners := make([]BamScanner, len(contamination))
	rejected := make([]bool, len(contamination))
	found := make([]bool, len(contamination))

	for c := 0; c < len(contamination); c++ {
		if err := contScanners[c].OpenBam(contamination[c]); err != nil {
			logger.Fatal(err)
		}
		reads_found[c] = 0
		reads_filtered[c] = 0
	}

	header, err := ReadBamHeader(args.Sample)
	if err != nil {
		logger.Fatal(err)
	}

	out := BamWriter{}
	outfp, err := out.Open(args.Output)
	if err != nil {
		logger.Fatal(err)
	}

	io.WriteString(outfp, header)

	reads_kept := 0
	read_mates_kept := 0
	total_reads := 0
	total_read_mates := 0
	ercc := 0
	considered := 0
	too_short := 0
	too_diverged := 0

	err = func() error {
		defer scanner.Done()
		defer benchmark(startedAt, "processing")

		for {
			if total_reads > 0 && total_reads%100000 == 0 {
				kept_percent = float64(reads_kept) / float64(considered) * 100
				logger.Printf("considered %d out of %d so far, kept %0.1f%%\n", considered, total_reads, kept_percent)
			}
			if args.Limit > 0 && args.Limit == total_reads {
				return nil
			}

			// Set up flags for outcomes wrt each potential source of contamination.
			for c, _ := range contamination {
				rejected[c] = false
				found[c] = false
			}

			// Read the first mate in a paired end run.
			mate1, err := scanner.Record()
			if err != nil {
				return fmt.Errorf("failed to read from sample BAM: %v after %d lines", err, total_reads)
			}
			if scanner.Closed {
				return nil
			}
			scanner.Ratchet()
			read := mate1[0]
			total_reads++
			total_read_mates++

			// See if we have the second mate of this pair.
			mate2, err := scanner.Find(read)
			if err != nil {
				return fmt.Errorf("failed to read from sample BAM: %v after %d lines", err, total_reads)
			}
			if mate2 != nil {
				scanner.Ratchet()
				total_read_mates++
			}

			var mate1_len int
			var mate1_edit_dist int
			var mate2_len int
			var mate2_edit_dist int

			mate1_len, mate1_edit_dist, err = extract(mate1)
			if err != nil {
				return err
			}
			if args.Verbose {
				logger.Println("found read", read, "mate 1:")
				logger.Println(strings.Join(mate1, "\t"))
			}
			if mate2 != nil {
				mate2_len, mate2_edit_dist, err = extract(mate2)
				if err != nil {
					return err
				}
				if args.Verbose {
					logger.Println("found read", read, "mate 2:")
					logger.Println(strings.Join(mate2, "\t"))
				}
			}

			// Filter for ERCC if either mate is mapped to ERCC.
			if MatchesErcc(mate1, mate2) {
				ercc++
				if args.Verbose {
					logger.Println("ERCC, rejecting")
				}
				continue
			}

			if mate1_len < args.MinLength {
				// If we don't have mate2 or if it's also too short, we mark this pair as too short.
				if mate2 == nil || mate2_len < args.MinLength {
					if args.Verbose {
						logger.Println("too short, rejecting")
					}
					too_short++
					continue
				}
				if args.Verbose {
					logger.Println("promoting mate 2")
				}
				// Mate2 is okay, so we promote it to mate1, and forget mate2
				mate1_len = mate2_len
				mate1_edit_dist = mate2_edit_dist
				mate1 = mate2
				mate2 = nil
			}
			if mate2 != nil && mate2_len < args.MinLength {
				// We have a mate2, but it doesn't meet the min length criteria, just forget it.
				mate2 = nil
				if args.Verbose {
					logger.Println("mate 2 too short, forgetting")
				}
			}
			// We treate the filter for edit distance the same way as length.
			if mate1_edit_dist > args.MaxDist {
				if mate2 == nil || mate2_edit_dist > args.MaxDist {
					too_diverged++
					if args.Verbose {
						logger.Println("too divergent, rejecting")
					}
					continue
				}
				if args.Verbose {
					logger.Println("promothing mate 2")
				}
				// Mate2 is okay, so we promote it to mate1, and forget mate2
				mate1_len = mate2_len
				mate1_edit_dist = mate2_edit_dist
				mate1 = mate2
				mate2 = nil
			}
			if mate2 != nil && mate2_edit_dist > args.MaxDist {
				// We have a mate2, but it doesn't meet the max edit distance criteria, just forget it.
				mate2 = nil
				if args.Verbose {
					logger.Println("mate 2, too diverged, forgetting")
				}
			}

			// If we get this far it means the read met the preliminary filtering criteria.
			considered++

			// Compare agains the best score for the read pair.
			mate1_score := float64(mate1_len) - float64(mate1_edit_dist)*args.Penalty
			var mate2_score float64
			best_score := mate1_score
			best_len := mate1_len
			best_edit_dist := mate1_edit_dist

			if mate2 != nil {
				mate2_score = float64(mate2_len) - float64(mate2_edit_dist)*args.Penalty
				if mate2_score > mate1_score {
					best_score = mate2_score
					best_len = mate2_len
					best_edit_dist = mate2_edit_dist
					if args.Verbose {
						logger.Printf("mate 2 has better score (%f) than mate 1 (%f)\n", mate2_score, mate1_score)
					}
				}
			}

			// Reads in the sample BAM will be rejected if either mate in any of the
			// contamination BAM files maps better than in the sampel BAM file.
			was_rejected := false
			for c := 0; c < len(contamination); c++ {
				m := 0
				for {
					mate, err := contScanners[c].Find(read)
					if err != nil {
						logger.Fatal(err)
					}
					if mate == nil {
						// No more alignments for this read in this contamination mapping
						break
					}
					m++
					if args.Verbose {
						logger.Printf("found mapping %d for %s in %s\n", m, mate[0], contamination[c])
						logger.Println(strings.Join(mate, "\t"))
					}
					if !found[c] {
						found[c] = true
						reads_found[c]++
					}
					length, edit_dist, err := extract(mate)
					if err != nil {
						logger.Fatalf("failed to read from %s: %v", contamination[c], err)
					}
					if length >= args.MinLength {
						score := float64(length) - float64(edit_dist)*args.Penalty
						if args.Verbose {
							logger.Printf("mapping meets length criteria and has score %f\n", score)
						}
						if best_score <= score+args.Margin {
							if args.Verbose {
								logger.Println("mapping has better score")
							}
							if !rejected[c] {
								reads_filtered[c]++
								rejected[c] = true
								was_rejected = true
								if args.Verbose {
									logger.Printf("read %s with length %d and edit distance %d was rejected "+
										"with score %0.1f because in %s it had a score of %0.1f with length "+
										"%d and edit distance %d\n",
										read, best_len, best_edit_dist, best_score, contamination[c],
										score, length, edit_dist)
								}
							}
						} else {
							if args.Verbose {
								logger.Println("mapping has worse score")
							}
						}
					}
				}
			}
			if !was_rejected {
				// This read is okay, output it to the output BAM file.
				_, err := fmt.Fprintf(outfp, "%s\n", strings.Join(mate1, "\t"))
				if err != nil {
					return err
				}
				reads_kept++
				read_mates_kept++
				if mate2 != nil {
					_, err := fmt.Fprintf(outfp, "%s\n", strings.Join(mate2, "\t"))
					if err != nil {
						return err
					}
					read_mates_kept++
				}
				if args.Verbose {
					logger.Printf("kept read %s with length %d and edit distance %d and score %0.1f\n",
						read, best_len, best_edit_dist, best_score)
				}
			}
		}
	}()
	if err != nil {
		logger.Fatal(err)
	}

	outfp.Close()
	out.Wait()

	logger.Println("Preliminary filtering:")
	if args.Ercc {
		erccPerc := float64(ercc) / float64(total_reads) * 100
		logger.Printf("filtered out %d ERCC reads (%0.1f%%) before comparing to contamination\n", ercc, erccPerc)
	}

	shortPerc := float64(too_short) / float64(total_reads) * 100
	logger.Printf("filtered out %d reads (%0.1f%%) becase their alignment was too short\n", too_short, shortPerc)
	divergedPerc := float64(too_diverged) / float64(total_reads) * 100
	logger.Printf("filtered out %d reads (%0.1f%%) becase they were too diverged\n", too_diverged, divergedPerc)

	logger.Printf("%d reads remaining after preliminary filtering\n", considered)
	logger.Println("Contamination filtering:")
	for c, cont := range contamination {
		n := reads_filtered[c]
		perc := float64(n) / float64(considered) * 100
		found_perc := float64(reads_found[c]) / float64(considered) * 100
		logger.Printf("found %d of %d reads in %s (%0.1f%%)\n", reads_found[c], considered, cont, found_perc)
		logger.Printf("rejected %d of %d reads from %s (%0.1f%%)\n", reads_filtered[c], considered, cont, perc)
	}

	kept_percent = float64(reads_kept) / float64(considered) * 100
	total_percent := float64(reads_kept) / float64(total_reads) * 100
	logger.Printf("kept %d of %d reads (%0.1f%%), which is %0.1f%% of the %d reads that met preliminary filtering\n",
		reads_kept, total_reads, total_percent, kept_percent, considered)
	total_mates_percent := float64(read_mates_kept) / float64(total_read_mates) * 100
	logger.Printf("kept %d of %d read mates (%0.1f%%)", read_mates_kept, total_read_mates, total_mates_percent)
	input_mates_per_pair := float64(total_read_mates) / float64(total_reads)
	output_mates_per_pair := float64(read_mates_kept) / float64(reads_kept)
	logger.Printf("observed %0.1f mates/read on the input end and %0.1f mates/read on the output end\n",
		input_mates_per_pair, output_mates_per_pair)

	logger.Println("machine parsable stats:")
	stats := []int{
		total_reads,
		total_read_mates,
		ercc,
		too_short,
		too_diverged,
		considered,
		reads_kept,
		read_mates_kept,
	}
	stats = append(stats, reads_found...)
	stats = append(stats, reads_filtered...)
	statsStr := "stats"
	for _, s := range stats {
		statsStr += fmt.Sprintf("\t%d", s)
	}
	logger.Println(statsStr)
}
