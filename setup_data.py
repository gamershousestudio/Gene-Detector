"""

This file exists to save genomes in order to analyze them in main.py.
Please run file and paste genome info from chosen database as per instructions.

"""

import sqlite3
import os

def write_genome(path, bases):
    # Connects to database
    connection = sqlite3.connect(path)
    cursor = connection.cursor()

    # Creates a table for genome
    cursor.execute("CREATE TABLE IF NOT EXISTS genome(base_sequence TEXT)")

    # Inserts bases
    cursor.execute("INSERT OR REPLACE INTO genome VALUES (?)", (bases,))

    # Passes data
    connection.commit()
    connection.close()

def get_genes(text):
    # Loops through all features and looks for keyword "gene"
    # Following code relies on patterns used in text files by genome.jp
    forward_genes = []

    for i in range(len(text)):
        if text[i:i+4] == "gene":
            if text[i+4].isdigit():
                # Finds the number for the start
                start = 0
                start_index = 0

                for j in range(i+4, len(text)):
                    if text[j] == ".": # Signifies end of start
                        start = int(text[i+4:j])
                        start_index = j
                        break

                stop = 0

                for j in range(start_index+2, len(text), 1):
                    if text[j] == "/": # Signifies end of start
                        stop = text[start_index+2:j]
                        break

                forward_genes.append((start, stop))

    return forward_genes

    # TODO: ADD PROTEIN SCANNING FOR REVERSED GENES

def write_genes(path, forward_genes, backward_genes):

    # Forward strand proteins
    with sqlite3.connect(path) as connection:

        cursor = connection.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS forward_proteins(starts INTEGER, stops INTEGER)")

        for start, stop in forward_genes:
            cursor.execute("INSERT INTO forward_proteins VALUES (?, ?)", (start, stop))

    # Temp debug
    connection = sqlite3.connect(path)
    cursor = connection.cursor()
    cursor.execute("SELECT starts, stops FROM forward_proteins")

def create_genome():
    # Accesses dataset folder
    folder = "Genomes"
    os.makedirs(folder, exist_ok=True)  # <-- creates folder if it doesn't exist

    # Makes a .db for the organism to map
    name = input("What is the name of the organism:\n").replace(' ', '_')
    db_path = os.path.join(folder, f"{name}-genome.db")

    print(f'\nPlease paste full genome document')

    # Loops through all lines on the pasted dataset
    lines = []
    while True:
        line = input()
        if "//" in line:
            break
        lines.append(line)

    # Joins all lines together
    bases = ''.join(lines)

    gene_text = "" # First section contains proteins (at least on genome.jp)

    # Checks if entire document pasted, if so looks for start
    if "ORIGIN" in bases:
        gene_text = bases.split("ORIGIN")[0].split("FEATURES")[1].replace(" ", "").replace("\n", "")
        bases = bases.split("ORIGIN")[1]

    # Removes all invalid characters
    bases = ''.join(c.upper() for c in bases if c.lower() in {'a', 't', 'c', 'g', 'n'})

    write_genome(db_path, bases)

    # Protein mapping(if applicable)
    if gene_text != '':
        forward = get_genes(gene_text)
        write_genes(db_path, forward, [])

if __name__ == '__main__':
    create_genome()