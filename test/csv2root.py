#!/usr/bin/env python3

import csv
import tqdm
import copy
import argparse
import os
from pyjetty.mputils import perror, pinfo, pwarning, treewriter
import ROOT
ROOT.gROOT.SetBatch(True)

def make_dict(scols, row):
	d = {}
	for i, s in enumerate(scols):
		d[s] = float(row[i])
	return d

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('input', type=str, default='')
	parser.add_argument('--nev', type=int, default=-1)
	args = parser.parse_args()
	finame = args.input
	f = open(finame)
	csv_f = csv.reader(f)

	rfout = ROOT.TFile(finame.replace('.csv', '.root'), 'recreate')
	rfout.cd()
	tw = treewriter.RTreeWriter(tree_name='emissions', fout=rfout)

	pinfo('working on a file', finame)
	irow = 0
	scols = []
	nmax = args.nev
	for row in tqdm.tqdm(csv_f):
		if irow < 1:
			scols = copy.deepcopy(row)
			irow = irow + 1
			continue
		else:
			tw.fill_branch('e', make_dict(scols, row))
			tw.fill_tree()
			if nmax > 0 and irow >= nmax:
				break
			irow = irow + 1
	rfout.Write()
	pinfo('file written:', rfout.GetName())
	pinfo('done.')


if __name__ == '__main__':
	main()