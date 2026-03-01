# Stuff to import
import math
import os
import sqlite3

import matplotlib.pyplot as plt
import numpy as np

# Genome information
bacterial = True

codon_table = { # Copied from chatgpt
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
    stops = {"TAA", "TAG", "TGA"}  # Which codons correspond to stop codons
    starts = {"ATG", "ATA", "GTG"} # Main start codons

    potential_genes = []

    # Seperates into three reading frames
    for frame in range(3):
        # Loops through each position in each reading frame
        for pos in range(frame, len(bases) - 2, 3):
            codon = ''.join(bases[pos:pos + 3])

            # Checks for start codon
            if codon in starts:
                # Looks for stop codon
                for stop in range(pos + 3, len(bases) - 2, 3):
                    stop_codon = ''.join(bases[stop:stop + 3])

                    if stop_codon in stops:
                        potential_genes.append((pos, stop + 3))
                        break

    return potential_genes

def remove_nested(potential_genes, confidence):
    # Sorts gene list by length
    potential_genes.sort(key=lambda x: x[1] - x[0], reverse=True)

    filtered = []
    filtered_confidence = []

    # Looks through all the start and stops
    for i, (start, stop) in enumerate(potential_genes):
        contained = False  # Marks it as not contained by default
        for fstart, fstop in filtered:

            # Flags it if it is nested inside another gene
            if start >= fstart and stop <= fstop:
                contained = True
                break

        # Filters out any not nested genes
        if not contained:
            filtered.append((start, stop))
            filtered_confidence.append(confidence[i])
        elif confidence[i] > .5:
            filtered.append((start, stop))
            filtered_confidence.append(confidence[i])

    potential_genes = filtered

    return potential_genes, confidence

def shine_dalgarno(potential_genes, bases, confidence_factor):
    confidence = [0 for i in range(len(potential_genes))]

    # Partials
    partials = {"GAGG", "GGAGGT", "AGGA", "GGAG", "GAGGT", "TAAGG"}

    if bacterial:
        for i, (start, stop) in enumerate(potential_genes):
            upstream = ''.join(bases[max(0, start - 10):start])

            if "AGGAGG" in upstream:
                confidence[i] += confidence_factor
            elif True in [partial in upstream for partial in partials]: # Checks for partials
                confidence[i] += confidence_factor

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

    # Checks for genes where 3rd reading frame has most gc content
    for i, (start, stop) in enumerate(potential_genes):
        gc1 = 0
        gc2 = 0
        gc3 = 0

        for codon in range(start, stop, 3):
            if bases[codon] == "G" or bases[codon] == "C":
                gc1 += 1
            if bases[codon+1] == "G" or bases[codon+1] == "C":
                gc2 += 1
            if bases[codon+2] == "G" or bases[codon+2] == "C":
                gc3 += 1

        if gc1 < gc3 and gc2 < gc3:
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

    # Increases confidence on genes with more preffered codons and vise versa
    for i, (start, stop) in enumerate(potential_genes):
        score = 0
        codon_count = 0

        for j in range(start, stop, 3):
            codon = bases[j:j + 3]

            if codon in codon_freq:
                score += codon_freq[codon]
                codon_count += 1

        if codon_count > 0:
            confidence[i] += score / codon_count

    return confidence

def kozak(potential_genes, bases, confidence_factor_positive, confidence_factor_negative):
    confidence = [0 for i in range(len(potential_genes))]

    # (gcc)gccRccATGG
    for i, (start, stop) in enumerate(potential_genes):
        if start < 10: continue # Not likely; but just to be safe

        inst_seq = ''.join(bases[start-9:start+1])

        if inst_seq[-1] == "G":
            confidence[i] += confidence_factor_positive
        else:
            confidence[i] -= confidence_factor_negative

        if inst_seq[6] == "G" or inst_seq[6] == "A":
            confidence[i] += confidence_factor_positive
        else:
            confidence[i] -= confidence_factor_negative

        if inst_seq[4:] == "GCCACCATGG":
            confidence[i] += confidence_factor_positive*2 # Exact kozak sequence

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

    # Subtract score from genes ending with least common stop codon
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

def check_genes(potential_genes, bases, bacterial, threshold):
    shine_dalgarno_confidence = shine_dalgarno(potential_genes, bases, 0.3) if bacterial else [0 for i in
                                                                                               range(len(potential_genes))]
    gc_confidence = gc_comparison(potential_genes, bases, .3, .1, .48)
    codon_bias_confidence = codon_bias_check(potential_genes, bases, .3, .1)
    kozak_confidence = kozak(potential_genes, bases, .3, .1)
    stop_distribution_confidence = stop_distribution(potential_genes, bases, .1)

    confidence = [(shine_dalgarno_confidence[i] + gc_confidence[i] + codon_bias_confidence[i] + kozak_confidence[i] + stop_distribution_confidence[i])
                  for i in range(len(potential_genes))]

    potential_genes, confidence = remove_nested(potential_genes, confidence)

    # Filters out super small genes
    filtered_genes = []
    filtered_confidence = []

    for i, (start, stop) in enumerate(potential_genes):
        if stop - start < 20:
            pass
        else:
            if stop - start < len(bases) / 1000: # If all genes failed this requernment only ~.1% of the genome would be genes; 10x smaller than the average percent in eukaryotes
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

    return potential_genes, confidence

# Graphing

def graph_line(length, genome, chromosome):
    plt.figure()
    plt.hlines(y=0, xmin=0, xmax=length, color="black")
    plt.xlim(0, length)

    plt.yticks([])
    plt.xlabel("Base Position")
    plt.title(f"{genome} chromosome {chromosome}")

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

if __name__ == "__main__":

    # Primary loop for accessing data
    while True:
        # Genome choice
        options = os.listdir("Genomes")

        options = [option for option in options if ".db" in option]
        options = [option.replace("-genome.db", "") for option in options]
        options = [option.replace("_", " ") for option in options]

        print("Please choose a genome to access:")
        for i, option in enumerate(options):
            print(f"{option} [{i + 1}]")

        try:
            genome = options[int(input("Enter index of chosen genome: ")) - 1].replace(" ", "_") + "-genome.db"
        except ValueError or IndexError or KeyError:
            print("Please enter valid index.")
            continue

        path = os.path.join("Genomes", genome)

        connection = sqlite3.connect(path)
        cursor = connection.cursor()

        # Gets table in dataset
        cursor.execute("SELECT chromosome_number FROM genome")
        tables = cursor.fetchall()

        # Chromosome choice
        access = input(f"Enter chromosome index to access(1 - {len(tables)}): ") if len(tables) > 1 else 1

        # Checks for valid choice
        passed = True
        try:
            access = int(access)
        except ValueError:
            passed = False

        if (access > len(tables) or access < 1) and passed:
            print("Invalid index.")

        # Checks passed; gives info on chromosome
        else:
            # Gets chromosome base sequence
            chromosome = access
            bases = cursor.execute("SELECT base_sequence FROM genome WHERE chromosome_number = ?", (chromosome,)).fetchall()[0][0]
            length = len(bases)

            print(f"chromosome: {access}")
            print(f"bases: {length}")

            graph_line(length, genome.replace("_", " ").removesuffix("-genome.db"), access)

            # Gets reverse complement for bases
            complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
            reverse_complement = "".join([complement[b] for b in reversed(bases) if b is not "N"])

            potential_genes = get_genes(bases)
            reverse_potential_genes = get_genes(reverse_complement)

            # Gene checks
            threshold = 1.2

            potential_genes, confidence = check_genes(potential_genes, bases, bacterial, threshold)
            reverse_potential_genes, reverse_confidence = check_genes(reverse_potential_genes, reverse_complement, bacterial, threshold)

            # Makes reversed genes start/stop accurate
            reverse_potential_genes = [(length - stop, length - start)for start, stop in reverse_potential_genes]

            # Plots gene map
            graph_genes(potential_genes, confidence, 1)
            graph_genes(reverse_potential_genes, reverse_confidence, -1)

                # Loops through all codons to assign amino acids
                # amino_acids = []

                # for j in range(int(start/3), int(stop/3)):
                    # amino_acids.append(codon_table[codons[int(j)]])


                # print(f"total bps in gene {i}: {len(amino_acids*3)}")
                # print(f"amino acids in gene {i}: {amino_acids}\n")

                # Updates current letter
                # i += 1

            # Plots TATA boxes
            #for i in range(len(bases)-3):
            #    if bases[i:i+6] == "TATAAA":
            #        plt.scatter(i, 0, color="green", s=10)

            plt.show()

            print(f"{len(potential_genes) + len(reverse_potential_genes)} genes mapped.")
            print(f"Genes(start, stop): forward {potential_genes} reverse {reverse_potential_genes}")

            if input("Press ENTER to exit, or any key then enter to view other chromosome\n") == '':
                break