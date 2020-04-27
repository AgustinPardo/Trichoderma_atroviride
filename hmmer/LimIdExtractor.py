#!/usr/bin/env python3

from Bio import SearchIO
import numpy as np
import argparse

def getID(string):
    try: return string.split('|')[1]
    except IndexError: return string

def filePrinter(limFile, IDfile, searchFile, threshold):
    lims=open(limFile, 'w')
    IDs=open(IDfile, 'w')
    for seq in searchFile:
        saveID=False
        for domain in seq:
            if domain.bitscore>=threshold:
                saveID=True
                print(getID(seq.id), domain.env_start, domain.env_end, file=lims)
        if saveID: print(getID(seq.id), file=IDs)
    lims.close()
    IDs.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", nargs='+', help="hmmsearch output file")
    parser.add_argument("-th", nargs='+', help="domain bitscore threshold")
    parser.add_argument("-lim", nargs='+', help="name for the domain limits file")
    parser.add_argument("-id", nargs='+', help="name for the ID file")
    args=parser.parse_args()
    search=SearchIO.read(args.i[0], 'hmmer3-text')
    filePrinter(args.lim[0], args.id[0], search, int(args.th[0]))
