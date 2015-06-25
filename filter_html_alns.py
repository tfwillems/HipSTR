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
        self.counts = {}
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
            self.counts[self.coord] = collections.defaultdict(int)
        self.counts[self.coord][data.upper()] += 1

class FilteredHTMLOutputter(HTMLParser):
    def __init__(self, skip_columns, left_trim):
        self.skip_columns = skip_columns
        self.left_trim    = left_trim
        HTMLParser.__init__(self)
        
    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.coord = 0

        output       = "<" + tag 
        coord_change = 1
        for i in xrange(len(attrs)):
            if attrs[i][0] == "colspan":
                output += " %s=\"%d\""%(attrs[i][0], int(attrs[i][1])-self.left_trim)
            else:
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
    nlines = 0
    for line in sys.stdin:
        lines.append(line)
        if not line.startswith("<tr") or "samplename" in line:
            continue
        parser.feed(line.strip())
        nlines += 1
    bad_cols = set()
    for key, vals in parser.counts.items():
        if len(vals) == 1 and "*" in vals:
            bad_cols.add(key)

    if nlines != 1:
        fixed_points = map(lambda x: x[0], sorted(filter(lambda x: len(x[1]) == 1 and (x[1].values()[0] == 1), parser.counts.items())))
        if len(fixed_points) >= 1 and fixed_points[0] == 1:
            trim_index = 1
            while trim_index < len(fixed_points) and fixed_points[trim_index] == fixed_points[trim_index-1]+1:
                trim_index += 1
            trim_index = max(0, trim_index-10)
            for i in xrange(0, trim_index):
                bad_cols.add(fixed_points[i])
            left_trim = 0 if trim_index == 0 else fixed_points[trim_index-1]
        else:
            left_trim = 0

        if len(fixed_points) >= 1 and fixed_points[-1] == sorted(parser.counts.keys())[-1]:
            trim_index = len(fixed_points)-2
            while trim_index >= 0 and fixed_points[trim_index]+1 == fixed_points[trim_index+1]:
                trim_index -= 1
            trim_index = min(len(fixed_points), trim_index+10)
            for i in xrange(trim_index, len(fixed_points)):
                bad_cols.add(fixed_points[i])
    else:
        left_trim = 0

    # Output each line after removing columns only containing a *
    output = FilteredHTMLOutputter(bad_cols, left_trim)
    for line in lines:
        if not line.startswith("<tr") or "samplename" in line:
            sys.stdout.write(line)
            continue
        output.feed(line.strip())
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
