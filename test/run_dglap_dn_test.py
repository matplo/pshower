#!/usr/bin/env python3

import argparse
import os
import dn
import tqdm

# long int RANDOMSEED = atoi(argv[1]);
# long int NUMEVENTS = atoi(argv[2]);
# double Q = atof(argv[3]);
# double Rmax = atof(argv[4]);
# long int Flavor = atoi(argv[5]);
# long int ProgramMode = atoi(argv[6]);

NUMEVENTSNoteWorking=1000
NUMEVENTSSave=50000
InvertQgluon=-2.0
InvertQquark=-0.5
SetNumInverts=10000000
NumTries=10000
TheNumEbins=800

import csv

def dict_from_MapSD(e):
	d = {}
	for _k in e:
		d[_k] = e[_k] 
	return d

def write_emissions_to_a_file(fname, emissions):
	if len(emissions) > 0:
		fieldnames = [k for k in emissions[0]]
	else:
		return
	if not os.path.exists(fname):
		with open(fname, 'w', newline='') as csvfile:
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
			writer.writeheader()

	dn.dump_emissions_csv(fname)
	# with open(fname, 'a', newline='') as csvfile:
	# 	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
	# 	for e in emissions:
	# 		d = dict_from_MapSD(e)
	# 		writer.writerow(d)

def main():
	# 101 ${nev} 250 1.5708 1 1
	parser = argparse.ArgumentParser(description='dglap DN implemenation', prog=os.path.basename(__file__))
	parser.add_argument('--initfile', type=str, default="Initialize_Parton_Shower.txt")
	parser.add_argument('--seed', type=int, default=10101)
	parser.add_argument('--nev', type=int, default=NUMEVENTSSave)
	parser.add_argument('-Q', type=float, default=250)
	parser.add_argument('--Rmax', type=float, default=1.5708)
	parser.add_argument('--flavor', type=int, default=1)
	parser.add_argument('--gg-only', help='force ->gg split', action='store_true', default=False)
	parser.add_argument('--program-mode', type=int, default=1)
	parser.add_argument('--filenameLeadingJetSpectra', default='', type=str)
	parser.add_argument('--filenameEventWideSpectra', default='', type=str)
	parser.add_argument('--filenameEventWideSpectraLogBin', default='', type=str)
	parser.add_argument('--csv', help='csv output for emissions', default='', type=str)

	args = parser.parse_args()

	if args.gg_only:
		dn.ExtraParameters.instance().set_flag("ggonly")
	csv_output_file = None
	if len(args.csv) > 0:
		if os.path.exists(args.csv):
			os.remove(args.csv)
		dn.ExtraParameters.instance().set_flag("get_emissions")

	params = dn.doubleArray(10)
	NumRadii = dn.longintp()
	NumRadii.assign(1)
	JetRadii = dn.doubleArray(100)

	dn.LoadShowerParams(args.initfile, params, NumRadii, JetRadii);

    # We note the vales in params:
    # *(params) = Minimum Energy Fraction or Zcut;
    # *(params+1) = CA;
    # *(params+2) = NF;
    # *(params+3) = CF;
    # *(params+4) = Jet Energy;
    # *(params+5) = Jet Initial Opening Angle;
    # *(params+6) = Minimum Mass Scale;

	params_text = [
		'Minimum Energy Fraction or Zcut',
		'CA',
		'NF',
		'CF',
		'Jet Energy',
		'Jet Initial Opening Angle',
		'Minimum Mass Scale']

	params[4] = args.Q
	params[5] = args.Rmax

	print('[i] dumping params')
	for i in range(len(params_text)):
		print(' - ', params_text[i], '=', params[i])
	print('[i] NumRadii = ', NumRadii.value())

	if args.filenameLeadingJetSpectra == '':
		args.filenameLeadingJetSpectra = "pyLeadingJetSpectra_rseed{}_Q{}_Rmax{}_ProgramMode{}_Flavor{}_Angle.txt".format(args.seed, args.Q, args.Rmax, args.program_mode, args.flavor)
	if args.filenameEventWideSpectra == '':
		args.filenameEventWideSpectra = "pyInclusiveSpectra_rseed{}_Q{}_Rmax{}_ProgramMode{}_Flavor{}_Angle.txt".format(args.seed, args.Q, args.Rmax, args.program_mode, args.flavor)
	if args.filenameEventWideSpectraLogBin == '':
		args.filenameEventWideSpectraLogBin = "pyInclusiveSpectraLogBin_rseed{}_Q{}_Rmax{}_ProgramMode{}_Flavor{}_Angle.txt".format(args.seed, args.Q, args.Rmax, args.program_mode, args.flavor)

	print(args)
	# for i in range(NumRadii.value()):
	# 	print(JetRadii[i])

	emissions = dn.PSEmissionsList()
	DaughterEmissions = dn.PSEmissionsList()
	# not used here
	# EmissionsWithinJet = dn.PSEmissionsList()

	dn.PSemissions_init(emissions)
	dn.PSemissions_init(DaughterEmissions)
	# dn.PSemissions_init(EmissionsWithinJet)

	t = dn.doublep()
	AveT = dn.doublep()
	RND = dn.doublep()
	# not used 
	# Z = dn.doublep()
	CurrentWTAaxis = dn.doubleArray(2)

	dn.ReInitialize(emissions, DaughterEmissions, CurrentWTAaxis, args.flavor, t);
	AveT.assign(0)

	# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
	# //Build Splitting Function Inversion
	print("Computing Inversion of Splitting Functions.")
	NumInverts = SetNumInverts
	InvertPg = dn.doubleArray(NumInverts+1)
	InvertPq = dn.doubleArray(NumInverts+1)

	r = dn.reset_random_number_generator(args.seed)
	dn.BuildSplitFunctionInversion(InvertPg, InvertPq, NumInverts, params, r);

	LeadingJetSpectra = dn.doubleParray(NumRadii.value(), TheNumEbins)
	EventWideSpectra = dn.doubleParray(NumRadii.value(), TheNumEbins)
	EventWideSpectraLogBin = dn.doubleParray(NumRadii.value(), TheNumEbins)

	for i in tqdm.tqdm(range(args.nev)):
		t.assign(0)
		dn.JetFragmentation(emissions, DaughterEmissions, CurrentWTAaxis, JetRadii,
							NumRadii.value(), params, TheNumEbins, LeadingJetSpectra, EventWideSpectra,
							EventWideSpectraLogBin, i+1, t, InvertPg, InvertPq, NumInverts, args.program_mode, r);
		AveT.assign(AveT.value() + t.value());
		# //This wipes out all the EmissionsList's, reseeds emissions with a parton pointed at the z axis with energy fraction 1, and flavor specified
		dn.ReInitialize(emissions, DaughterEmissions, CurrentWTAaxis, args.flavor, t);

		_save_intermediate = False
		if (i+1) % NUMEVENTSSave == 0:
			_save_intermediate = True

		if args.nev < NUMEVENTSSave and (i+1 == args.nev):
			_save_intermediate = True

		if i < 1:
			continue

		if _save_intermediate:
			dn.WriteToDiskHistgramsFRAG(args.filenameLeadingJetSpectra, (i+1), NumRadii.value(),
				TheNumEbins, LeadingJetSpectra, JetRadii, dn.MinZBook,
				args.flavor, params, 0);
			dn.WriteToDiskHistgramsFRAG(args.filenameEventWideSpectra, (i+1), NumRadii.value(),
				TheNumEbins, EventWideSpectra, JetRadii, dn.MinZBook,
				args.flavor, params, 0);
			dn.WriteToDiskHistgramsFRAG(args.filenameEventWideSpectraLogBin, (i+1), NumRadii.value(),
				TheNumEbins, EventWideSpectraLogBin, JetRadii, dn.MinZBook,
				args.flavor, params, 0);
	
		if i > 0 and i % 1000 == 0:
			if len(args.csv) > 0:
				write_emissions_to_a_file(args.csv, dn.ExtraOutput.instance().emissions)
				dn.ExtraOutput.instance().emissions.clear()

	if len(args.csv) > 0:
		write_emissions_to_a_file(args.csv, dn.ExtraOutput.instance().emissions)
		dn.ExtraOutput.instance().emissions.clear()
		print('[i] csv file:', args.csv)
	print('[i] done.')


if __name__ == '__main__':
	main()

