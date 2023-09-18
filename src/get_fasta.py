#!/usr/bin/env python
#
# Copyright 2008 Gabriele Sales <gbrsales@gmail.com>

from __future__ import with_statement
from optparse import OptionParser
from sys import stdin, stdout
from vfork.fasta.reader import MultipleBlockStreamingReader
from vfork.io.util import safe_rstrip
from vfork.util import format_usage, exit, ignore_broken_pipe
import re

def identity_filter(m):
	return m

def id_filter(m):
	sep_rx = re.compile(r'[ \t]')
	def f(value):
		match = sep_rx.search(value)
		if match:
			value = value[:match.start()]
		return m(value)
	return f

class ExactMatch(object):
	__slots__ = ['text']
	
	def __init__(self, text):
		self.text = text
	
	def __call__(self, value):
		return self.text == value

class RegularExpressionMatch(object):
	__slots__ = ['rx']
	
	def __init__(self, rx):
		self.rx = re.compile(rx)
	
	def __call__(self, value):
		return self.rx.search(value) is not None

def WordRegularExpressionMatch(rx):
	return RegularExpressionMatch('(^|\W)%s(\W|$)' % rx)

class SingleLabelSearch(object):
	def __init__(self, match, reuse_match, writer):
		self.match = match
		self.reuse_match = reuse_match
		self.writer = writer
	
	def scan(self, reader):
		found = False
		for header, lines in reader:
			if self.match(header):
				self.writer.write(header, lines)
				found = True
				if not self.reuse_match:
					break

		return found

class MultipleLabelSearch(object):
	def __init__(self, matches, reuse_matches, writer):
		self.matches = matches
		self.reuse_matches = reuse_matches
		self.writer = writer
	
	def scan(self, reader):
		found = [ False ] * len(self.matches)
		for header, lines in reader:
			for idx, match in enumerate(self.matches):
				if match(header):
					self.writer.write(header, lines)
					if self.reuse_matches:
						found[idx] = True
					else:
						del self.matches[idx]
						del found[idx]
					break
			
			if len(self.matches) == 0:
				break
		
		return all(found)

class IndexSearch(object):
	def __init__(self, idx, writer):
		self.idx = idx
		self.writer = writer
	
	def scan(self, reader):
		pos = 0
		for header, lines in reader:
			if pos == self.idx:
				self.writer.write(header, lines)
				break
			else:
				pos += 1
		else:
			return False
		
		return True

class RangeSearch(object):
	def __init__(self, start, stop, writer):
		self.start = start
		self.stop = stop
		self.writer = writer
	
	def scan(self, reader):
		pos = -1
		for header, lines in reader:
			pos += 1
			if self.stop is not None and pos >= self.stop:
				break
			elif self.start is None or pos >= self.start:
				self.writer.write(header, lines)
		else:
			if self.start is not None and pos < self.start:
				return False
			elif self.stop is not None and pos < self.stop-1:
				return False
			else:
				return True
		
		return True

class FastaWriter(object):
	def __init__(self, print_headers):
		self.print_headers = print_headers
	
	def write(self, header, lines):
		if self.print_headers:
			stdout.write('>%s\n' % header)
		for line in lines:
			stdout.write(line + '\n')

def is_range(value):
	return len(value) and value[0] == '@'

def parse_range(value):
	value = value[1:]
	idx = value.find(':')
	if idx == -1:
		try:
			idx = int(value)
			if idx < 0: raise ValueError
		except ValueError:
			exit('Invalid index: %s' % value)
		
		return ('idx', idx)
	
	else:
		start = value[:idx]
		if len(start) == 0:
			start = None
		else:
			try:
				start = int(start)
				if start < 0: raise ValueError
			except ValueError:
				exit('Invalid start index: %s' % start)
		
		stop = value[idx+1:]
		if len(stop) == 0:
			stop = None
		else:
			try:
				stop = int(value[idx+1:])
				if stop <= start: raise ValueError
			except ValueError:
				exit('Invalid stop index: %s' % value[idx+1:])
		
		return ('range', start, stop)

def load_labels(filename):
	try:
		labels = set()
		with file(filename, 'r') as fd:
			for line in fd:
				label = safe_rstrip(line)
				if label in labels:
					exit('Duplicated label \'%s\' in file %s' % (label, filename))
				else:
					labels.add(label)
		
		if len(labels) == 0:
			exit('No label in file %s' % filename)

		return labels
	
	except IOError:
		exit('Error reading from file %s' % filename)

def main():
	parser = OptionParser(usage=format_usage('''
		%prog [OPTIONS] SPEC <FASTA
		
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
	'''))
	parser.add_option('-i', '--id-only', dest='id_only', action='store_true', default=False, help='match the header only up to the first space or tab')
	parser.add_option('-r', '--regexp', dest='regexp', action='store_true', default=False, help='match headers against a regular expression')
	parser.add_option('-w', '--word-regexp', dest='word_regexp', action='store_true', default=False, help='scan headers for a word matching the given regular expression')
	parser.add_option('-f', '--from-file', dest='from_file', action='store_true', default=False, help='read labels from PATH', metavar='PATH')
	parser.add_option('-n', '--no-headers', dest='no_headers', action='store_true', default=False, help='do not print headers')
	parser.add_option('-m', '--ignore-missing', dest='ignore_missing', action='store_true', default=False, help='do not signal an error if some of the requested blocks could not be found')
	parser.add_option('-p', '--read-all-input-file', dest='read_all_input_file', action='store_true', default=False, help='when all requested block are printed do not close the stdin pipe but continue reading in order to avoid pustream SIGPIPE')
	options, args = parser.parse_args()
	
	if len(args) != 1:
		exit('Unexpected argument number.')
	elif options.regexp and options.word_regexp:
		exit('Cannot use --regexp and --word-regexp together.')
	elif is_range(args[0]) and (options.regexp or options.word_regexp):
		exit('Cannot use --regexp or --word-regexp with indexes.')
	
	r = MultipleBlockStreamingReader(stdin, join_lines=False)
	w = FastaWriter(not options.no_headers)
	
	if is_range(args[0]):
		spec = parse_range(args[0])
		if spec[0] == 'idx':
			s = IndexSearch(spec[1], w)
		else:
			s = RangeSearch(spec[1], spec[2], w)
	
	else:
		if options.regexp:
			mc = RegularExpressionMatch
			reuse_match = True
		elif options.word_regexp:
			mc = WordRegularExpressionMatch
			reuse_match = True
		else:
			mc = ExactMatch
			reuse_match = False
		
		if options.id_only:
			fc = id_filter
		else:
			fc = identity_filter
		
		if options.from_file:
			s = MultipleLabelSearch([ fc(mc(l)) for l in load_labels(args[0]) ], reuse_match, w)
		else:
			s = SingleLabelSearch(fc(mc(args[0])), reuse_match, w)
	
	if not s.scan(r) and not options.ignore_missing:
		exit('Some FASTA blocks are missing.')
 
	
	if options.read_all_input_file:
		for line in stdin:
			pass

if __name__ == '__main__':
	ignore_broken_pipe(main)
