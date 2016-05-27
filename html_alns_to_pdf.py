# Simple class that processes HTML files created by HipSTR
# and removes any alignment positions in which all bases correspond
# to an insertion
# Useful when filtering HTML files for a subset of samples and insertions
# are no longer present

import collections
import sys
from HTMLParser import HTMLParser
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

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
        if self.coord not in self.counts:
            self.counts[self.coord] = collections.defaultdict(int)
        self.counts[self.coord][data.upper()] += 1


cell_width   = 8
cell_height  = 16
cell_padding = 1
text_colors  = {"locustd":"black", "Atd":"purple",   "Ctd":"blue",     "Gtd":"green",     "Ttd":"orange", "vtd":"gray", "-td":"red",
                "Areftd":"white",  "Creftd":"white", "Greftd":"white", "Treftd":"white",  "vreftd":"white",
                "areftd":"white",  "creftd":"white", "greftd":"white", "treftd":"white",  "snptd":"black", "inserttd": "black"}

fill_colors  = {"locustd":"N/A",   "Atd":"white",    "Ctd":"white",    "Gtd":"white",     "Ttd":"white", "vtd":"white", "-td":"white",
                "Areftd":"purple", "Creftd":"blue",  "Greftd":"green", "Treftd":"orange", "vreftd":"gray",
                "areftd":"purple", "creftd":"blue",  "greftd":"green", "treftd":"orange", "snptd": "#ffd700", "inserttd": "red"}

class FilteredPDFOutputter(HTMLParser):
    def __init__(self, skip_columns, left_trim, num_lines, num_columns, output_file):
        self.output_stream = open(output_file, "w")
        self.skip_columns = skip_columns
        self.left_trim    = left_trim
        HTMLParser.__init__(self)
        self.row_number   = 1
        self.output_stream.write("<svg width=\"%d\" height=\"%d\">\n"%((num_columns+1)*(cell_width+cell_padding), (num_lines+10)*(cell_height+cell_padding)))

    def handle_caption(self, line):
        self.row_number += 2
        text = line.strip().split("<caption>")[1].split("</caption>")[0].replace("\t", " ")
        x = 5
        y = (self.row_number-1)*(cell_height+cell_padding) + cell_padding
        res  = "<text x=\"%d\" y=\"%d\" font-family=\"%s\" font-size=\"%s\" font-weight=\"bold\" fill=\"%s\">%s</text>\n"%(x, y, "Courier", "15px", "black", text)
        self.output_stream.write(res)

    def handle_sample(self, line):        
        self.row_number += 2
        text = line.strip().split("</font>")[0]
        text = text.split("<font color=\"red\">")[1]
        x = 5
        y = (self.row_number-1)*(cell_height+cell_padding) + cell_padding
        res  = "<text x=\"%d\" y=\"%d\" font-family=\"%s\" font-size=\"%s\" fill=\"%s\">%s</text>\n"%(x, y, "Courier", "13px", "red", text)
        self.output_stream.write(res)
        
    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.coord       = 0
            self.row_number += 1
            self.viz_coord   = 0

        if tag != "td":
            return

        output       = "<g>"
        coord_change = 1
        num_cells    = 1
        self.fill_color = "black"
        self.text_color = "black"
        self.bold_font  = False
        self.spanner    = False
        self.missing_color = False
        for i in xrange(len(attrs)):
            if attrs[i][0] == "colspan":
                lt        = len(filter(lambda x: x < int(attrs[i][1]), self.skip_columns))
                num_cells = int(attrs[i][1]) - self.left_trim - lt
                coord_change = int(attrs[i][1])
                self.spanner = True
            elif attrs[i][0] == "class":
                if "td" in attrs[i][1]:
                    self.fill_color = fill_colors[attrs[i][1]]
                    self.text_color = text_colors[attrs[i][1]]
                    if attrs[i][1] == "locustd":
                        self.bold_font = True

                    if attrs[i][1] == "snptd" or attrs[i][1] == "inserttd":
                        self.missing_color = True

            elif attrs[i][0] == "style":
                self.fill_color = attrs[i][1].replace("background-color:","")
                
        x         = (self.viz_coord)*(cell_width+cell_padding) + cell_padding
        y         = (self.row_number-1)*(cell_height+cell_padding) + cell_padding
        width     = num_cells*cell_width + (num_cells-1)*cell_padding
        height    = cell_height
        if not self.spanner and self.fill_color != "white":
            output += "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\"/>"%(x, y, width, height, self.fill_color)

        self.coord_change = coord_change
        if (self.coord + coord_change) not in self.skip_columns:
            self.output_stream.write(output)
            self.skip = False
        else:
            self.skip = True

    def handle_endtag(self, tag):
        if tag == "td":
            for i in xrange(1, self.coord_change+1):
                if (self.coord + i) not in self.skip_columns:
                    self.viz_coord += 1

            self.coord += self.coord_change
            if not self.skip:
                self.output_stream.write("</g>\n")

    def handle_data(self, data):
        if not self.skip:
            x    = (self.viz_coord)*(cell_width+cell_padding) + cell_padding + cell_width/2.0
            y    = (self.row_number-1)*(cell_height+cell_padding) + cell_padding + 0.75*cell_height
            if self.missing_color:
                if data + "td" in text_colors:
                    self.text_color = text_colors[data+"td"]
            res  = "<text x=\"%d\" y=\"%d\" font-family=\"%s\" font-size=\"%s\" fill=\"%s\" text-anchor=\"middle\""%(x, y, "Courier", "13px", self.text_color)
            res += "" if not self.bold_font else " font-weight=\"800\""
            res += " alignment-baseline=\"central\">%s</text>"%(data)
            self.output_stream.write(res)


    def finish(self):
        self.output_stream.write("</svg>\n")
        self.output_stream.close()


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
        fixed_points = map(lambda x: x[0], sorted(filter(lambda x: len(x[1]) == 1 and (x[1].values()[0] == 1), parser.counts.items())))
        if len(fixed_points) >= 1 and fixed_points[0] == 1:
            trim_index = 1
            while trim_index < len(fixed_points) and fixed_points[trim_index] == fixed_points[trim_index-1]+1:
                trim_index += 1
            trim_index = max(0, trim_index-10)
            for i in xrange(0, trim_index):
                if parser.counts[trim_index][0] != "*":
                    bad_cols.add(fixed_points[i])
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
    num_columns = max(parser.counts.keys()) - len(bad_cols)
    output_path_prefix = sys.argv[1]
    output = FilteredPDFOutputter(bad_cols, left_trim, nlines, num_columns, output_path_prefix+".svg")
    for line in lines:
        if line.startswith("<style"):
            continue
        if "caption" in line:
            output.handle_caption(line)
            continue
        if "samplename" in line:
            output.handle_sample(line)
            continue
        
        if not line.startswith("<tr") or "samplename" in line:
            continue
        output.feed(line.strip())
        #sys.stdout.write("\n")
    output.finish()

    drawing  = svg2rlg(output_path_prefix + ".svg")
    renderPDF.drawToFile(drawing, output_path_prefix + ".pdf")


if __name__ == "__main__":
    main()
