#!/usr/bin/env python3

from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-lim", nargs='+', help="domain limits file")
parser.add_argument("-fas", nargs='+', help="fasta file with full sequence")
parser.add_argument("-o", nargs='+', help="fasta with domains")
args=parser.parse_args()


fasta=SeqIO.parse(args.fas[0], 'fasta')

matrixID=[]
matrixLim=[]

with open (args.lim[0], 'r') as inf:
    for item in inf:
        matrixLim.append(item.split()[1:])
        matrixID.append(item.split()[0])

def intMatrix(matrix):
    matrix2=[]
    for x in matrix:
        matrix2.append([int(x[0]), int(x[1])])
    return matrix2

matrixLim=intMatrix(matrixLim)

def getID(string):
    try: return string.split('|')[1]
    except IndexError: return string

with open (args.o[0], 'w') as outf:
    for record in fasta:
        for k, acc in enumerate(matrixID):
            if getID(record.id)==acc:
                print('>'+getID(record.id)+'/'+str(matrixLim[k][0])+'-'+str(matrixLim[k][1])+'\n'+record.seq[matrixLim[k][0]:matrixLim[k][1]], file=outf)
