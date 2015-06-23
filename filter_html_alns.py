#!/usr/local/bin/python

# Simple class that processes HTML files created by HipSTR
# and removes any alignment positions in which all bases correspond
# to an insertion
# Useful when filtering HTML files for a subset of samples and insertions
# are no longer present

import collections
import sys
from HTMLParser import HTMLParser

class HTMLCharCounter(HTMLParser):
    def __init__(self):
        self.counts = collections.defaultdict(int)
        HTMLParser.__init__(self)
        
    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.coord = 0

        coord_change = 1
        for i in xrange(len(attrs)):
            if tag == "td":
                if attrs[i][0] == "colspan":
                    coord_change = int(attrs[i][1])
        if tag == "td":
            self.coord += coord_change

    def handle_data(self, data):
        if data == ' ':
            return
        if self.coord not in self.counts:
            self.counts[self.coord] = set()
        self.counts[self.coord].add(data.upper())

class FilteredHTMLOutputter(HTMLParser):
    def __init__(self, skip_columns):
        self.skip_columns = skip_columns
        HTMLParser.__init__(self)
        
    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.coord = 0

        output       = "<" + tag 
        coord_change = 1
        for i in xrange(len(attrs)):
            output += " %s=\"%s\""%(attrs[i][0], attrs[i][1])
            if tag == "td":
                if attrs[i][0] == "colspan":
                    coord_change = int(attrs[i][1])
        if tag == "td":
            self.coord += coord_change
        output += ">"

        if tag != "td" or self.coord not in self.skip_columns:
            sys.stdout.write(output)
            self.skip = False
        else:
            self.skip = True

    def handle_endtag(self, tag):
        if not self.skip:
            sys.stdout.write("</" + tag + ">")

    def handle_data(self, data):
        if not self.skip:
            sys.stdout.write(data)


def main():
    lines = []

    # Read all lines in first pass and determine which columns
    # only have a * character (insertion)
    parser = HTMLCharCounter()
    for line in sys.stdin:
        lines.append(line)
        if not line.startswith("<tr") or "samplename" in line:
            continue
        parser.feed(line.strip())
    bad_cols = set()
    for key, vals in parser.counts.items():
        if len(vals) == 1 and "*" in vals:
            bad_cols.add(key)

    # Output each line after removing columns only containing a *
    output = FilteredHTMLOutputter(bad_cols)
    for line in lines:
        if not line.startswith("<tr") or "samplename" in line:
            sys.stdout.write(line)
            continue
        output.feed(line.strip())
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
