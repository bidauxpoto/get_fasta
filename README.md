# get_fasta

Usage: get_fasta [OPTIONS] SPEC <FASTA

Reads a FASTA file and retrieves the blocks corresponding to SPEC.

SPEC can be one of the following:
* LABEL - search for a block with the given label (see also -i, -r and
          -w modifiers);
* @IDX - extract the block at index IDX (zero-based).
* @START:STOP - extract the blocks at indexes between START and STOP
         (start position is inclusive, stop is exclusive).
         Indexes are zero-base and can be both omitted: a missing START
         implies 0; a missing STOP corresponds to the end of the input;
* FILENAME - used in conjuction with -f modifier. Labels are read from
         an external file, one per line.

Options:
```
  -h, --help            show this help message and exit
  -i, --id-only         match the header only up to the first space or tab
  -r, --regexp          match headers against a regular expression
  -w, --word-regexp     scan headers for a word matching the given regular
                        expression
  -f, --from-file       read labels from PATH
  -n, --no-headers      do not print headers
  -m, --ignore-missing  do not signal an error if some of the requested blocks
                        could not be found
  -p, --read-all-input-file
                        when all requested block are printed do not close the
                        stdin pipe but continue reading in order to avoid
                        pustream SIGPIPE
```
