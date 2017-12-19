# Simple class that processes HTML files created by HipSTR
# and removes any alignment positions in which all bases correspond
# to an insertion
# Useful when filtering HTML files for a subset of samples and insertions
# are no longer present

import collections
import sys

# Import HTMLParser, whose package name changed between python2 and python3
if sys.version_info[0] == 2:
    from HTMLParser import HTMLParser
elif sys.version_info[0] == 3:
    from html.parser import HTMLParser
else:
    exit("Unsupported python version %d"%(sys.version_info[0]))

class HTMLCharCounter(HTMLParser):
    def __init__(self):
        self.counts = {}
        HTMLParser.__init__(self)

    def handle_data(self, data):
        offset, bases = data.split(" ")
        offset = int(offset)
        for i in range(len(bases)):
            if offset not in self.counts:
                self.counts[offset] = collections.defaultdict(int)
            self.counts[offset][bases[i].upper()] += 1
            offset += 1

class FilteredHTMLOutputter(HTMLParser):
    def __init__(self, skip_columns, left_trim):
        self.skip_columns = skip_columns
        self.left_trim    = left_trim

        self.true_base_dict = {'H': 'A', 'I': 'C', 'J':'G', 'K':'T', 'L':'N'}
        self.color_dict     = {
            'A': "purple",
            'C': "blue",
            'G': "green",
            'T': "orange",
            'N': "purple",
            '-': "red"}
        HTMLParser.__init__(self)
        
    def handle_starttag(self, tag, attrs):
        if tag != "tr":
            exit("ERROR: Improperly formatted input file")
        sys.stdout.write("<tr>")
        self.reference = False
        for t,v in attrs:
            if t == "class" and v == "reference":
                self.reference = True

    def handle_endtag(self, tag):
        if tag != "tr":
            exit("ERROR: Improperly formatted input file")
        sys.stdout.write("</tr>")

    def handle_data(self, data):
        offset, bases = data.split(" ")
        offset        = int(offset)
        lt            = len(list(filter(lambda x: x < offset, self.skip_columns)))
        new_offset    = offset - self.left_trim - lt
        if new_offset > 0:
            sys.stdout.write("<td colspan=%d> </td>"%(new_offset))
        coord = offset

        for i in range(len(bases)):
            if self.reference:
                if bases[i] == "*":
                    output="<td class=\"%sreftd\">%s</td>"%("v", bases[i])
                else:
                    output="<td class=\"%sreftd\">%s</td>"%(bases[i], bases[i])
            else:
                if bases[i] == 'x':
                    output="<td class=\"%std\">%s</td>"%("spacer", bases[i])
                elif bases[i] == '*':
                    output="<td class=\"%std\">%s</td>"%("v", bases[i])
                elif bases[i] in ['A', 'C', 'G', 'T', 'N', '-']:
                    output="<td class=\"%std\">%s</td>"%(bases[i], bases[i])
                elif bases[i] in ['a', 'c', 'g', 't', 'n']:
                    true_base = bases[i].upper()
                    output="<td class=\"inserttd\"><font color=\"%s\">%s</font></td>"%(self.color_dict[true_base], true_base)
                elif bases[i] in ['H', 'I', 'J', 'K', 'L']:
                    true_base = self.true_base_dict[bases[i]]
                    output="<td class=\"snptd\"><font color=\"%s\">%s</font></td>"%(self.color_dict[true_base], true_base)
                else:
                    print(bases[i])
                    exit("Invalid base character in file: " + bases[i])

            if coord not in self.skip_columns:
                sys.stdout.write(output)
            coord += 1


def write_header():
    header="""
    <style type='text/css'>

     .ref {
      color: white;
      font-family: Courier;
     }

     td {
      text-align:center;
      vertical-align:middle;
      font-family: Courier;
      font-size: 13px;
     }

     .locustd {
      font-style: italic;
      color: black;
     }

     .snptd {
       background-color: gold;
     }

     .inserttd {
      background-color: red;
     }

     .spacertd {
      color: white;
     }

     .reftable {
      color: white;
     }

     .readtable {
      font-weight: normal;
     }

     caption {
      background: #dbb768;
      color:black;
      font-weight: bold;
      font-size: 1.1em;
      text-align: left;
     }

     .Atd {
      color: purple;
     }

     .Ctd {
      color: blue;
     }

     .Gtd {
      color: green;
     }

     .Ttd {
      color: orange;
     }

     .vtd {
      color: gray;
     }

     .-td {
      color: red;
     }

     .Areftd {
      background-color: purple;
     }

     .Creftd {
      background-color: blue;
     }

     .Greftd {
      background-color: green;
     }

     .Treftd {
      background-color: orange;
     }

     .vreftd {
      background-color: gray;
     }

     </style>
    """
    print(header.replace("\n", " "))

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
        elif len(vals) == 2 and "*" in vals and " " in vals:
            bad_cols.add(key)

    if nlines != 1:
        fixed_points = list(map(lambda x: x[0], sorted(filter(lambda x: len(x[1]) == 1 and (list(x[1].values())[0] == 1), list(parser.counts.items())))))
        if len(fixed_points) >= 1 and fixed_points[0] == 1:
            trim_index = 1
            while trim_index < len(fixed_points) and fixed_points[trim_index] == fixed_points[trim_index-1]+1:
                trim_index += 1
            trim_index = max(0, trim_index-10)
            for i in range(0, trim_index):
                if parser.counts[trim_index][0] != "*":
                    bad_cols.add(fixed_points[i])
        left_trim = 0

        if len(fixed_points) >= 1 and fixed_points[-1] == sorted(parser.counts.keys())[-1]:
            trim_index = len(fixed_points)-2
            while trim_index >= 0 and fixed_points[trim_index]+1 == fixed_points[trim_index+1]:
                trim_index -= 1
            trim_index = min(len(fixed_points), trim_index+10)
            for i in range(trim_index, len(fixed_points)):
                bad_cols.add(fixed_points[i])
    else:
        left_trim = 0

    # Output each line after removing columns only containing a *
    write_header()
    output = FilteredHTMLOutputter(bad_cols, left_trim)
    for line in lines:
        if not line.startswith("<tr") or "samplename" in line:
            sys.stdout.write(line)
            continue
        output.feed(line.strip())
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
