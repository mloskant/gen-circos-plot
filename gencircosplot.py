from argparse import ArgumentParser
import pybedtools
import pycircos
import RNA
import collections
import cv2
import os

# Command line flags
parser = ArgumentParser()
parser.add_argument("-g", "--genes", dest="genes", help="The genes in BED format", required=True)
parser.add_argument("-b", "--bindings", dest="bindings", nargs='+',
                    help="The binding sites of proteins in BED format",
                    required=False)
parser.add_argument("-fi", "--fasta", dest="fasta", help="The chromosomal sequences as a fasta file", required=True)
parser.add_argument("-o", "--output", dest="output", help="Folder in which all pngs are saved", required=True,
                    default="myCircos")

args = parser.parse_args()
genes = str(args.genes)
fasta_file = args.fasta
folder = str(args.output)

os.mkdir(folder)

if args.bindings and len(args.bindings) > 5:
    parser.error("A maximum of 5 arguments may be specified for --bindings.")


def convert_dot_bracket(structure, gene_name):
    """
    Converts dot-bracket structure for circos

    Parameters
    ----------
    structure : str
        Dot-bracket structure of a sequence
    gene_name: str
        Name of the gene

    Returns
    -------
    start_list: list
        Start data location of linked data
    end_list: list
        End data location of linked data
    """
    startpos = 1
    currentpos = 1
    start_list = []
    end_list = []
    counter = 0
    while startpos < len(structure):
        while currentpos <= len(structure):
            # The positions with points are skipped
            if structure[currentpos - 1] == ".":
                currentpos += 1
                continue
            # If there is an open bracket at the current position, the counter is incremented
            if structure[currentpos - 1] == "(":
                counter += 1
                currentpos += 1
                continue
            # If there is a closed bracket at the current position, the counter is decremented
            if structure[currentpos - 1] == ")":
                counter -= 1
                currentpos += 1
            # If the counter is 0, a base pairing has been found.
            # More precisely, the matching closed bracket to the open bracket "()" has been found.
            if counter == 0:
                end_list = end_list + [(gene_name, currentpos - 1, currentpos - 1, 749)]
                start_list = start_list + [(gene_name, startpos, startpos, 749)]
                startpos += 1
                # The next starting position is being sought
                while startpos <= len(structure) and structure[startpos - 1] != "(":
                    startpos += 1
                currentpos = startpos
                break

    return start_list, end_list


def tighten_linked_data(start_list, end_list):
    """
    Tights linked data to avoid having individual strands in the circos plot that actually belong together

    Parameters
    ----------
    start_list : list
        Start data location of linked data
    end_list: list
        End data location of linked data

    Returns
    -------
    start_list_belt: list
        Start data location of linked data
    end_list_belt:
        End data location of linked data
    """
    start_list_belt = []
    end_list_belt = []
    current_position = 0
    # The complete lists consisting of the start positions/end positions of the linked data are run through
    while current_position < len(start_list) - 1:
        _, start1, _, raxis1_start = start_list[current_position]
        _, start2, _, _ = start_list[current_position + 1]
        _, end1, _, raxis1_end = end_list[current_position]
        _, end2, _, _ = end_list[current_position + 1]
        # If the next start position is one greater than the current start position
        # and the corresponding end positions also suit, then they are combined as a belt
        if start1 + 1 == start2 and end2 + 1 == end1:
            start1_belt = start1
            end2_belt = end1
            current_position += 1
            # The belt is extended until the positions are no longer consecutive,
            # then the new start and end positions of the belts are saved
            extension = 2
            if current_position + 1 < len(start_list):
                _, start2, _, _ = start_list[current_position + 1]
                _, end2, _, _ = end_list[current_position + 1]
                while start1 + extension == start2 and end2 + extension == end1 and current_position < len(
                        start_list) - 1:
                    current_position += 1
                    if current_position + 1 < len(start_list):
                        _, start2, _, _ = start_list[current_position + 1]
                        _, end2, _, _ = end_list[current_position + 1]
                    extension += 1
            _, start2_belt, _, _ = start_list[current_position]
            _, end1_belt, _, _ = end_list[current_position]
            start_list_belt = start_list_belt + [(gene_id, start1_belt, start2_belt, raxis1_start)]
            end_list_belt = end_list_belt + [(gene_id, end1_belt, end2_belt, raxis1_end)]
        # If the positions are not consecutive, the next position is considered
        else:
            current_position += 1
    return start_list_belt, end_list_belt


def find_bindings_on_gene(bindings, gene_name, gene_chrom, gene_chromstart, gene_chromend):
    """
    Finds and saves the protein bindings on which the gene is

    Parameters
    ----------
    bindings: str
        The binding sites of proteins in BED format
    gene_name: str
        Name of the gene
    gene_chrom: str
        Chromosome on which the gene is located
    gene_chromstart: int
        Start of the sequence of the gene
    gene_chromend: int
        End of the sequence of the gene

    Returns
    -------
    arcdata_dict: dict
        A dictionary of binding site positions and widths
    """
    arcdata_dict = collections.defaultdict(dict)
    bindings_bed = pybedtools.BedTool(bindings)
    # Filter out all irrelevant binding sites
    for feature in bindings_bed:
        if feature.chrom == gene_chrom:
            if feature.start <= gene_chromend and feature.end >= gene_chromstart:
                if feature.strand == gene_strand:
                    start = int(feature.start)
                    width = int(feature.end) - int(feature.start)
                    if gene_name not in arcdata_dict:
                        arcdata_dict[gene_name]["positions"] = []
                        arcdata_dict[gene_name]["widths"] = []
                    # Save position on the gene
                    arcdata_dict[gene_name]["positions"].append(start - gene_chromstart)
                    # Save length of the binding site
                    arcdata_dict[gene_name]["widths"].append(width)
    return arcdata_dict


def hex_to_rgb(hex_color):
    """
    Converts a hexadecimal color string to an RGB tuple.

    Parameters
    ----------
    hex_color: str
        A string representing a color in hexadecimal format.

    Returns
    -------
    rgb: tuple
        A tuple containing three integers that represent the red, green, and blue components of the color.
    """
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return rgb


# Load the BED file
bed_file = pybedtools.BedTool(genes)

# Run through the entries in the BED file
for feature in bed_file:
    gene_sequence = ""
    gene_chrom = feature.chrom  # chromosome on which the gene is located
    gene_chromstart = feature.start  # start coordinate on the chromosome
    gene_chromend = feature.end  # end coordinate on the chromosome
    gene_strand = feature.strand  # strand orientation
    gene_id = feature.name  # id
    # Extract the sequence from the FASTA file
    sequence = bed_file.sequence(fi=fasta_file, name=True, s=True)
    with open(sequence.seqfn, 'r') as file:
        for line in file:
            # save the sequence at the appropriate line until the next gene begins
            if line.startswith(">" + gene_id):
                try:
                    line = next(file)
                    gene_sequence = gene_sequence + line
                    line = next(file)
                    while not (line.startswith(">")):
                        gene_sequence = gene_sequence + line
                        line = next(file)
                except StopIteration:
                    break

    # Folding of the RNA sequence
    # Package viennarna version: 2.6.4, RNAfold
    sequence = gene_sequence.strip()
    structure, _ = RNA.fold(sequence)

    # Circos plot
    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle
    circle = Gcircle(figsize=(20, 20))
    # The circle itself
    label = "This is the circos of " + gene_id
    arc = Garc(arc_id=gene_id, size=len(sequence), raxis_range=(750, 750), labelposition=200, labelsize=30,
               label_visible=True, edgecolor="black", facecolor="black", linewidth=1.5, label=label)
    circle.add_garc(arc)
    # The arc rectangle with start and end angle of the circos plot
    circle.set_garcs(0, 363)

    # Markings on the circle, here in an interval of 100
    # BED file starts with 0
    if len(sequence) > 4000:
        ticklabels = list(range(0, len(sequence) + 1, 1000))
    elif len(sequence) == 2500:
        ticklabels = list(range(0, len(sequence), 100))
    else:
        ticklabels = list(range(0, len(sequence) + 1, 100))
    circle.tickplot(gene_id, raxis_range=(750, 760), tickwidth=2.0,
                    tickpositions=ticklabels, ticklabels=ticklabels,
                    ticklabeldirection="outer")

    # Chord plot for each base pairing
    # Converting dot-bracket structure for circos
    color = ["lightgreen", "limegreen", "mediumseagreen", "seagreen", "green", "darkgreen"]
    start_list, end_list = convert_dot_bracket(structure, gene_id)
    start_list_belt, end_list_belt = tighten_linked_data(start_list, end_list)
    color_variable = 0
    for i in range(0, len(start_list_belt)):
        circle.chord_plot(start_list_belt[i], end_list_belt[i], linewidth=0.2, facecolor=color[color_variable])
        color_variable += 1
        if color_variable > len(color) - 1:
            color_variable = 0

    # Save the plot as specified in the command line
    circle.figure.savefig(folder + "/" + gene_id + ".png")

    # Barplot for protein bindings
    if args.bindings is not None:
        # Position of the bindings
        raxis_range_start = 760
        raxis_range_end = 805

        bindings = args.bindings
        color = ['#1E90FF', '#FFC0CB', '#8a2be2', '#ff4040', '#ff8c00']
        s = 0
        for binding in bindings:
            arcdata_dict = find_bindings_on_gene(binding, gene_id, gene_chrom, gene_chromstart, gene_chromend)
            for key in arcdata_dict:
                circle.barplot(key, data=[1] * len(arcdata_dict[key]["positions"]),
                               positions=arcdata_dict[key]["positions"],
                               width=arcdata_dict[key]["widths"], raxis_range=[raxis_range_start, raxis_range_end],
                               edgecolor="black", facecolor=color[s], linewidth=0.5)
            s += 1
            raxis_range_start += 55
            raxis_range_end = raxis_range_start + 45
        # Save the plot as specified in the command line
        circle.figure.savefig(folder + "/" + gene_id + ".png")
        i = 0
        # Legende
        # Position of the legend/rectangle
        position_y = 1750
        position_y2 = 1770
        for binding in bindings:
            # Load the image
            image = cv2.imread(folder + "/" + gene_id + ".png")

            legend_text = binding
            legend_position = (1600, position_y2)

            # rectangle (BGR color)
            rgb = hex_to_rgb(color[i])
            rectangle_color = (rgb[2], rgb[1], rgb[0])
            rectangle_top_left = (1550, position_y)
            rectangle_bottom_right = (1590, position_y2)
            cv2.rectangle(image, rectangle_top_left, rectangle_bottom_right, rectangle_color, thickness=cv2.FILLED)

            # Font, size, colour and thickness for the text
            font = cv2.FONT_HERSHEY_SIMPLEX
            font_scale = 0.5
            font_color = (0, 0, 0)  # Schwarz
            font_thickness = 1
            cv2.putText(image, legend_text, legend_position, font, font_scale, font_color, font_thickness)

            # Save the image with the added legend
            cv2.imwrite(folder + "/" + gene_id + ".png", image)

            # position for next text and rectangle
            position_y += 35
            position_y2 += 35
            # next color
            i += 1
