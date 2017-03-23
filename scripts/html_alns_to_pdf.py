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
        self.max_coord = 0
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
        self.max_coord = max(self.max_coord, self.coord)


cell_width   = 8
cell_height  = 13
cell_padding = 2
text_colors  = {"locustd":"white", "Atd":"purple",   "Ctd":"blue",     "Gtd":"green",     "Ttd":"orange", "vtd":"gray", "-td":"red",
                "Areftd":"white",  "Creftd":"white", "Greftd":"white", "Treftd":"white",  "vreftd":"white",
                "areftd":"white",  "creftd":"white", "greftd":"white", "treftd":"white",  "snptd":"black", "inserttd": "black"}

fill_colors  = {"locustd":"N/A",   "Atd":"white",    "Ctd":"white",    "Gtd":"white",     "Ttd":"white", "vtd":"white", "-td":"white",
                "Areftd":"purple", "Creftd":"blue",  "Greftd":"green", "Treftd":"orange", "vreftd":"gray",
                "areftd":"purple", "creftd":"blue",  "greftd":"green", "treftd":"orange", "snptd": "#ffd700", "inserttd": "#ff4d4d"}

class FilteredPDFOutputter(HTMLParser):
    def __init__(self, num_lines, num_columns, output_file):
        self.output_stream = open(output_file, "w")
        HTMLParser.__init__(self)
        self.row_number   = 1
        self.output_stream.write("<svg width=\"%d\" height=\"%d\">\n"%((num_columns+1)*(cell_width+cell_padding), (num_lines+10)*(cell_height+cell_padding)))

    def handle_caption(self, line):
        self.row_number += 2
        text = line.strip().split("<caption>")[1].split("</caption>")[0].replace("\t", " ")
        x = 5
        y = (self.row_number-1)*(cell_height+cell_padding) + cell_padding
        res  = "<text x=\"%d\" y=\"%d\" font-family=\"%s\" font-size=\"%s\" font-weight=\"bold\" fill=\"%s\">%s</text>\n"%(x, y, "Courier", "15px", "white", text)
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
                num_cells = int(attrs[i][1])
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
            output += "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\"/>"%(x, y-cell_padding/2, width+cell_padding, height+cell_padding, self.fill_color)

        self.coord_change = coord_change
        self.output_stream.write(output)

    def handle_endtag(self, tag):
        if tag == "td":
            for i in xrange(1, self.coord_change+1):
                self.viz_coord += 1

            self.coord += self.coord_change
            self.output_stream.write("</g>\n")

    def handle_data(self, data):
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
    lines    = []
    max_cols = 0
    nlines   = 0
    parser   = HTMLCharCounter()
    for line in sys.stdin:
        lines.append(line)
        if not line.startswith("<tr") or "samplename" in line:
            nlines += ("samplename" in line)
            continue
        parser.feed(line.strip())
        nlines += 1

    num_columns = parser.max_coord + 5
    output_path_prefix = sys.argv[1]
    output = FilteredPDFOutputter(nlines+5, num_columns, output_path_prefix+".svg")
    for line in lines:
        if line.startswith("<style"):
            continue
        if "<caption" in line:
            output.handle_caption(line)
            continue
        if "samplename" in line:
            output.handle_sample(line)
            continue
        if not line.startswith("<tr") or "samplename" in line:
            continue
        output.feed(line.strip())
    output.finish()

    drawing  = svg2rlg(output_path_prefix + ".svg")
    renderPDF.drawToFile(drawing, output_path_prefix + ".pdf")


if __name__ == "__main__":
    main()
