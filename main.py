# Stuff to import
import math
import os
import sqlite3

import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np

# Genome information
threshold = 4
minimum_gene_length = 100

codon_table = { # Copied from ChatGPT
    "TTT":"Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu",
    "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",
    "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met",
    "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",
    "TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser",
    "CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "TAT":"Tyr", "TAC":"Tyr", "TAA":"*", "TAG":"*",
    "CAT":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "AAT":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "GAT":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "TGT":"Cys", "TGC":"Cys", "TGA":"*", "TGG":"Trp",
    "CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "AGT":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
    "GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"
}

# Checks

def get_genes(bases):
    stops = {"TAA", "TAG", "TGA"}
    starts = {"ATG", "ATA", "GTG", "TTG"}

    bases = bases.upper().replace("\n","").replace(" ","")  # clean input
    potential_genes = []

    for frame in range(3):
        pos = frame
        while pos < len(bases) - 2:
            codon = bases[pos:pos+3]
            if codon in starts:
                # search for stop in the same frame
                stop_pos = pos + 3
                while stop_pos < len(bases) - 2:
                    stop_codon = bases[stop_pos:stop_pos+3]
                    if stop_codon in stops:
                        potential_genes.append((pos, stop_pos+3))
                        break
                    stop_pos += 3
                pos += 3  # move to next codon in this frame
            else:
                pos += 3  # move to next codon in this frame

    return potential_genes

def shine_dalgarno(potential_genes, bases, confidence_factor):
    confidence = [0 for i in range(len(potential_genes))]

    # Partials
    partials = {"GAGG", "GGAGGT", "AGGA", "GGAG", "GAGGT", "TAAGG"}
    for i, (start, stop) in enumerate(potential_genes):
        upstream = ''.join(bases[max(0, start - 12):start-4])

        if "AGGAGG" in upstream:
            confidence[i] += confidence_factor
        elif True in [partial in upstream for partial in partials]: # Checks for partials
            confidence[i] += confidence_factor/4

    return confidence

def gc_comparison(potential_genes, bases, confidence_factor_positive, confidence_factor_negative, gc_goal):
    confidence = [0 for i in range(len(potential_genes))]

    # Checks for gc content
    for i, (start, stop) in enumerate(potential_genes):
        seq = bases[start:stop]
        gc = (seq.count("G") + seq.count("C")) / len(seq)

        # Genes generally have more GC contents
        if (gc_goal - .1) < gc < (gc_goal + .1):
            confidence[i] += confidence_factor_positive
        elif stop - start > 1000 and (gc_goal - .2) < gc < (gc_goal + .2): # More lenient on larger genes
            confidence[i] += confidence_factor_positive

    # Checks for genes where 3rd reading frame has most gc content
    for i, (start, stop) in enumerate(potential_genes):
        gc1 = 0
        gc2 = 0
        gc3 = 0

        c = 0

        for codon in range(start, stop, 3):
            if bases[codon] == "G" or bases[codon] == "C":
                gc1 += 1
            if bases[codon+1] == "G" or bases[codon+1] == "C":
                gc2 += 1
            if bases[codon+2] == "G" or bases[codon+2] == "C":
                gc3 += 1

            c += 1

        # Normalizes results
        for base in "ATCG":
            gc1 /= c
            gc2 /= c
            gc3 /= c

        periodicity_score = gc3 - (gc1 + gc2) / 2
        if periodicity_score > 0:
            confidence[i] += confidence_factor_positive
        else:
            confidence[i] -= confidence_factor_negative
    return confidence

def codon_bias_check(potential_genes, bases, confidence_factor_positive, confidence_factor_negative):
    confidence = [0 for i in range(len(potential_genes))]

    # Gets total codon usage in organism
    codon_freq = {}

    for start, stop in potential_genes:
        if stop - start > 1000:  # Only check longer genes(just to be safer)
            for j in range(start, stop, 3):
                codon = bases[j:j + 3]
                codon_freq[codon] = codon_freq.get(codon, 0) + 1

    # Normalizes the values
    total = sum(codon_freq.values())
    codon_freq = {k: v / total for k, v in codon_freq.items()}

    # Increases confidence on genes with more preferred codons and vise versa
    for i, (start, stop) in enumerate(potential_genes):
        score = 0
        codon_count = 0

        for j in range(start, stop, 3):
            codon = bases[j:j + 3]

            if codon in codon_freq:
                score += codon_freq[codon]
                codon_count += 1

        if codon_count > 0:
            confidence[i] += score * confidence_factor_positive / codon_count + (stop - start) / 1000
        else:
            confidence[i] -= confidence_factor_negative;

    return confidence

def stop_distribution(potential_genes, bases, confidence_factor_negative):
    confidence = [0 for i in range(len(potential_genes))]

    # Stop codon distribution
    taa = 0
    tag = 0
    tga = 0
    for start, stop in potential_genes:
        if bases[stop - 3:stop] == "TAA":
            taa += 1
        elif bases[stop - 3:stop] == "TAG":
            tag += 1
        elif bases[stop - 3:stop] == "TGA":
            tga += 1

    # Subtract score from genes ending with the least common stop codon
    if min(taa, tag, tga) == taa:
        for i, (start, stop) in enumerate(potential_genes):
            if bases[stop - 3:stop] == "TAA":
                confidence[i] -= confidence_factor_negative
    elif min(taa, tag, tga) == tga:
        for i, (start, stop) in enumerate(potential_genes):
            if bases[stop - 3:stop] == "TGA":
                confidence[i] -= confidence_factor_negative
    elif min(taa, tag, tga) == tag:
        for i, (start, stop) in enumerate(potential_genes):
            if bases[stop - 3:stop] == "TAG":
                confidence[i] -= confidence_factor_negative

    return confidence

def alternate_stops(potential_genes, bases, confidence_factor_positive, confidence_factor_negative):
    stops = {"TAA", "TAG", "TGA"}  # Which codons correspond to stop codons

    confidence = [0 for i in range(len(potential_genes))]

    for i, (start, stop) in enumerate(potential_genes):
        frame1 = 0
        frame2 = 0

        for j in range(start + 1, stop - 2, 3):  # frame +1
            if bases[j:j + 3] in stops:
                frame1 += 1

        for j in range(start + 2, stop - 2, 3):  # frame +2
            if bases[j:j + 3] in stops:
                frame2 += 1

        frame1 = frame1 / (stop - start)
        frame2 = frame2 / (stop - start)

        if frame1 > .02 or frame2 > .02:
            confidence[i] += confidence_factor_positive * (frame1 + frame2)
        else:
            confidence[i] -= confidence_factor_negative

    return confidence

# TODO: CODON ADAPTATION INDEX -- Wait --

# TODO: AMINO ACID ENTROPY CHECK -- Complete --
def amino_acid_entropy(potential_genes, bases, confidence_factor_positive, confidence_factor_negative, min_val, max_val):
    confidence = [0 for i in range(len(potential_genes))]

    for i, (start, stop) in enumerate(potential_genes):
        # Measures frequency of amino acids used
        codon_freq = {}

        t = 0 # Will be used to normalize values

        # Counts distribution of amino acids
        for codon_index in range(start, stop, 3):
            codon = bases[codon_index:codon_index + 3]
            amino_acid = codon_table[codon]

            if amino_acid in codon_freq:
                codon_freq[amino_acid] += 1
            elif not amino_acid == "*":
                codon_freq[amino_acid] = 1

            t += 1

        # Normalizes values
        for amino_acid in codon_freq:
            codon_freq[amino_acid] /= t

        # Computes shannon entropy
        # Entropy = - sum (p * log2(p))

        entropy = - sum([p * math.log2(p) for a, p in codon_freq.items()])

        # Scores based on randomness of ORF
        if min_val < entropy < max_val:
            confidence[i] += confidence_factor_positive
        elif stop - start < 1000:
            confidence[i] -= confidence_factor_negative * ((3 - entropy) if entropy < 3 else 2) / (stop - start)

    return confidence

# TODO: HELP PREVENT FALSE GENE DETECTION WITHIN GENES ON 'ATG' STARTS
def filter_stops(potential_genes, confidence, bases):
    best = {}

    # Organize genes by order
    pair = [(gene, confidence[i]) for i, gene in enumerate(potential_genes)]
    pair.sort(key = lambda x: x[0][0])

    potential_genes = [p[0] for p in pair]
    confidence = [p[1] for p in pair]

    # Loop through all genes
    for gene, conf in zip(potential_genes, confidence):
        start, stop = gene

        if stop not in best:
            best[stop] = (gene, conf)
        elif conf > best[stop][1]:
            if bases[start:start+3] == "ATG":
                if conf > 2*best[stop][1] and start > best[stop][0][0]: # Commonly exists in real genes; must be SUPER confidence
                    best[stop] = (gene, conf)
                elif start < best[stop][0][0]:
                    best[stop] = (gene, conf)
            else:
                best[stop] = (gene, conf)
    genes = [v[0] for v in best.values()]
    confidences = [v[1] for v in best.values()]

    return genes, confidences

# TODO: START CODON PREFERENCE
def start_codon_preference(potential_genes, bases, confidence_factor_positive, confidence_factor_negative):
    confidence = [0 for i in range(len(potential_genes))]

    for i, (start, stop) in enumerate(potential_genes):
        start_codon = bases[start:start+3]

        # Checks start codons
        if start_codon == "ATG": # Most common
            confidence[i] += confidence_factor_positive
        elif start_codon == "GTG": # Somewhat uncommon(~15%)
            confidence[i] -= confidence_factor_negative/2
        elif start_codon == "TTG": # Super uncommon
            confidence[i] -= confidence_factor_negative
        else: # All others very rare
            confidence[i] -= confidence_factor_negative * 2

    return confidence

# TODO: LENGTH-SCALED LIKELIHOOD (exponentially longer orfs are punished; less likely to exist)

def base_bias(potential_genes, bases, confidence_factor_positive, confidence_factor_negative, min_score):
    confidence = [0 for i in range(len(potential_genes))]

    for i, (start, stop) in enumerate(potential_genes):
        # Base usage counts
        frame1 = { "A" : 0, "T" : 0, "C" : 0, "G" : 0}
        frame2 = {"A" : 0, "T" : 0, "C" : 0, "G" : 0}
        frame3 = {"A" : 0, "T" : 0, "C" : 0, "G" : 0}

        # Counts distribution of bases
        for codon in range(start, stop, 3):
            frame1[bases[codon]] += 1
            frame2[bases[codon + 1]] += 1
            frame3[bases[codon + 2]] += 1

        # Normalizes all values
        for j in range(4):
            frame1[list(frame1.keys())[j]] /= ((stop - start) / 3)
            frame2[list(frame2.keys())[j]] /= ((stop - start) / 3)
            frame3[list(frame3.keys())[j]] /= ((stop - start) / 3)

        # Computes total difference in values between two frames
        frame1_score = (abs(frame1["A"] - frame2["A"]) +
                        abs(frame1["T"] - frame2["T"]) +
                        abs(frame1["C"] - frame2["C"]) +
                        abs(frame1["G"] - frame2["G"]))

        frame2_score = (abs(frame2["A"] - frame3["A"]) +
                        abs(frame2["T"] - frame3["T"]) +
                        abs(frame2["C"] - frame3["C"]) +
                        abs(frame2["G"] - frame3["G"]))

        frame3_score = (abs(frame3["A"] - frame1["A"]) +
                        abs(frame3["T"] - frame1["T"]) +
                        abs(frame3["C"] - frame1["C"]) +
                        abs(frame3["G"] - frame1["G"]))

        total = frame1_score + frame2_score + frame3_score

        if total > min_score:
            confidence[i] += confidence_factor_positive * total
        else:
            confidence[i] -= confidence_factor_negative


    return confidence

def check_genes(potential_genes, bases, threshold):
    shine_dalgarno_confidence = shine_dalgarno(potential_genes, bases, 1)
    gc_confidence = gc_comparison(potential_genes, bases, 1, 1, .5)
    codon_bias_confidence = codon_bias_check(potential_genes, bases, 20, 3)
    stop_distribution_confidence = stop_distribution(potential_genes, bases, .05)
    alternate_stops_confidence = alternate_stops(potential_genes, bases, 30, .2)
    base_bias_confidence = base_bias(potential_genes, bases, .3, .3, .5)
    entropy_confidence = amino_acid_entropy(potential_genes, bases, .3, .5, 2.5, 4.15)
    start_codon_confidence = start_codon_preference(potential_genes, bases, .3, .3)

    #print(f"gene {potential_genes[1737]} \ngc: {gc_confidence[1737]}, codon bias: {codon_bias_confidence[1737]}, kozak: {kozak_confidence[1737]}, stop distribution: {stop_distribution_confidence[1737]}")

    confidence = [(shine_dalgarno_confidence[i] + gc_confidence[i] + codon_bias_confidence[i] + stop_distribution_confidence[i]
                   + alternate_stops_confidence[i] + base_bias_confidence[i] + entropy_confidence[i] + start_codon_confidence[i])
                  for i in range(len(potential_genes))]

    # Used to track certain genes
    debug = True
    genes = [(711, 2636), (930, 2636)]

    for gene in genes:
        if debug:
            try:
                i = potential_genes.index((gene[0]-1, gene[1])) # Program offsets start by 1 by default
                print(f"Gene: {gene}")
                print(f"Confidence: {confidence[i]}")
                print(
                    f"Shine Dalgarno: {shine_dalgarno_confidence[i]}; GC distribution: {gc_confidence[i]}; codon bias: {codon_bias_confidence[i]}; stop usage: {stop_distribution_confidence[i]}; "
                    f"alternate stops: {alternate_stops_confidence[i]}; base bias: {base_bias_confidence[i]}; amino acid entropy: {entropy_confidence[i]}; start codon preference: {start_codon_confidence[i]}")
            except ValueError:
                print(f"Gene {gene} not found")

    # Filters out minimal genes
    filtered_genes = []
    filtered_confidence = []

    for i, (start, stop) in enumerate(potential_genes):
        if stop - start < minimum_gene_length:
            pass
        else:
            if stop - start < len(bases) / 1000: # If all genes failed this requirement only ~.1% of the genome would be genes; 10x smaller than the average percent in eukaryotes
                confidence[i] -= .2

            filtered_genes.append(potential_genes[i])
            filtered_confidence.append(confidence[i])

    potential_genes = filtered_genes
    confidence = filtered_confidence

    # Removes any genes with too low a confidence
    filtered_genes = []
    filtered_confidence = []

    for i, val in enumerate(confidence):
        if val < threshold:
            pass
        else:
            filtered_genes.append(potential_genes[i])
            filtered_confidence.append(confidence[i])

    potential_genes = filtered_genes
    confidence = filtered_confidence

    potential_genes, confidence = filter_stops(potential_genes, confidence, bases)

    return potential_genes, confidence

# Graphing

def graph_line(length, genome):
    plt.figure()
    plt.hlines(y=0, xmin=0, xmax=length, color="black")
    plt.xlim(0, length)

    plt.yticks([])
    plt.xlabel("Base Position")
    plt.title(genome)

    plt.scatter([0, 0, length], [-1, 1, 1], color="black", s=0)

def graph_genes(lines, alphas, y_scale):
    for i, (start, stop) in enumerate(lines):
        plt.scatter(start, 0, c="red", s=0)
        plt.scatter(stop, 0, c="blue", s=0)

        # Labels genes
        if start % 3 == 0:
            plt.hlines(y=.1*y_scale, xmin=start, xmax=stop, color="red", alpha=min(alphas[i] ** 2, 1))
        elif start % 3 == 1:
            plt.hlines(y=.2*y_scale, xmin=start, xmax=stop, color="blue", alpha=min(alphas[i] ** 2, 1))
        elif start % 3 == 2:
            plt.hlines(y=.3*y_scale, xmin=start, xmax=stop, color="brown", alpha=min(alphas[i] ** 2, 1))
        # plt.annotate(text=str(i),
        # xy=((start+stop)/2, 0),
        # xytext=((start+stop)/2, (.1 if i % 2 == 1 else -.1)),
        # ha="center")

        print(f"{i} of {len(lines)} genes mapped." + (
            "  Confidence = " + str(alphas[i]) if alphas[i] > 0 else ""))
        print(f"Gene mapped from"f" {start} to {stop}")

# Other

def quick_scan(): # Taken from setup_data.py
    print(f'\nPlease paste bases in chromosome {i + 1}.')
    print("Press ENTER on empty line when complete.\n")

    # Loops through all lines on the pasted dataset
    lines = []
    while True:
        line = input()
        if line == '':
            break
        lines.append(line)

    # Joins all lines together
    chromosome = ''.join(lines)

    # Removes all invalid characters
    chromosome = ''.join(c.upper() for c in chromosome if c.lower() in {'a', 't', 'c', 'g', 'n'})

    return chromosome

if __name__ == "__main__":

    passed = True

    # Primary loop for accessing data
    while True:
        # Genome choice
        options = os.listdir("Genomes")

        options = [option for option in options if ".db" in option]
        options = [option.replace("-genome.db", "") for option in options]
        options = [option.replace("_", " ") for option in options]

        print("Please choose a genome to access:")
        print("Quick Scan [0]")
        for i, option in enumerate(options):
            print(f"{option} [{i + 1}]")

        do_quick_scan = False

        bases = ""

        try:
            choice = int(input("Enter index of chosen genome: "))
            if choice != 0: genome = options[choice - 1].replace(" ", "_") + "-genome.db"
            else:
                bases = quick_scan()
                do_quick_scan = True

                try:
                    threshold = float(input(f"Enter accuracy threshold: "))
                except ValueError:
                    pass


        except ValueError or IndexError or KeyError:
            print("Please enter valid index.")
            continue

        if not do_quick_scan:
            path = os.path.join("Genomes", genome)

            connection = sqlite3.connect(path)
            cursor = connection.cursor()

            # Gets table in dataset
            cursor.execute("SELECT chromosome_number FROM genome")
            tables = cursor.fetchall()

            # Chromosome choice
            access = input(f"Enter chromosome index to access(1 - {len(tables)}): ") if len(tables) > 1 else 1
            try:
                threshold = float(input(f"Enter accuracy threshold: "))
            except ValueError:
                pass

            # Checks for valid choice
            passed = True
            try:
                access = int(access)
            except ValueError:
                passed = False

            if (access > len(tables) or access < 1) and passed:
                print("Invalid index.")

            # Checks passed; gives info on chromosome
        if passed:
            if not do_quick_scan:
                # Gets chromosome base sequence
                chromosome = access
                bases = cursor.execute("SELECT base_sequence FROM genome WHERE chromosome_number = ?",
                                       (chromosome,)).fetchall()[0][0]

            length = len(bases)
            print(f"bases: {length}")

            if not do_quick_scan:
                graph_line(length, ''.join([genome.replace("_", " ").removesuffix("-genome.db") + (f" Chromosome {access}" if len(tables) > 1 else "")]))
            else:
                graph_line(length, "Quick Scan")

            # Gets reverse complement for bases
            complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
            reverse_complement = "".join([complement[b] if b != "N" else "N" for b in reversed(bases)])

            potential_genes = get_genes(bases)
            reverse_potential_genes = get_genes(reverse_complement)

            potential_genes, confidence = check_genes(potential_genes, bases, threshold)
            reverse_potential_genes, reverse_confidence = check_genes(reverse_potential_genes, reverse_complement, threshold)

            # Makes reversed genes start/stop accurate
            reverse_potential_genes = [(length - stop, length - start)for start, stop in reverse_potential_genes]

            # Plots gene map
            graph_genes(potential_genes, confidence, 1)
            graph_genes(reverse_potential_genes, reverse_confidence, -1)

            plt.show(block=False)
            plt.pause(0.001) # Time for rendering to occur

            potential_genes.sort(key=lambda x: x[0])

            print(f"{len(potential_genes) + len(reverse_potential_genes)} genes mapped.")
            print(f"Genes(start, stop): forward {[(start+1, stop) for start, stop in potential_genes]} reverse {[(start+1, stop) for start, stop in reverse_potential_genes]}")

            if input("Press ENTER to exit, or any key then enter to view other organisms\n") == "":
                break