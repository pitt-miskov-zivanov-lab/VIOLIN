#Created by Difei May 2023

# Usage
###############
# python reach_to_biorecipes.py JSON_FILE
###############

# Function
###############
# process a JSON file using INDRA engine, and return outputs in BioRECIPEs format
###############

import argparse
import sys,os
import numpy as np
from indra.belief import BeliefEngine
from indra.statements import statements
from indra.statements import Agent
from indra.tools import assemble_corpus as ac
import requests
from indra.sources import indra_db_rest as idr
from indra.sources import reach
from indra import config
from indra.tools import assemble_corpus as ac
from indra.belief import SimpleScorer
import time
import json
import pandas as pd

# TODO: organize these functions and develop a class to handle the translation

BioRECIPE_reading_col = ['Regulator Name', 'Regulator Type', 'Regulator Subtype',
						 'Regulator HGNC ID', 'Regulator Database', 'Regulator ID',
						 'Regulator Compartment', 'Regulator Compartment ID',
						 'Regulated Name', 'Regulated Type', 'Regulated Subtype',
						 'Regulated HGNC ID', 'Regulated Database', 'Regulated ID',
						 'Regulated Compartment', 'Regulated Compartment ID',
						 'Sign', 'Connection Type', 'Mechanism', 'Site',
						 'Cell Line', 'Cell Type', 'Tissue Type', 'Organism',
						 'Score', 'Source', 'Statements', 'Paper IDs']

def get_indra_stmts_by_reach(file_name, citation=None, organism_priority=None):
	"""Return indra statements by processing the given REACH json file

	Parameters
	----------
	NOTE: required parameters by indra interface
	file_name : str
        The name of the json file to be processed.
    citation : Optional[str]
        A PubMed ID passed to be used in the evidence for the extracted INDRA
        Statements. Default: None
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
	"""

	rp = reach.process_json_file(file_name)
	return rp

def indra_stmts_reading_output(processor, output_file, indra_stats=False, flute=False, violin=False):
	"""This function outputs the reading file from indra statements

	Parameters
	----------
	processor:
		Different types of processor containing the extracted INDRA Statements
	output_file: xlsx
		The name of reading file containing the extracted interactions
	indra_stats : boolean
		Lists of indra statements in spreadsheet
	flute : boolean
		Lists of interactions in FLUTE supported format
	violin : boolean
		Reading output supported by VIOLIN old version

	Returns
	-------
	By default, it will output the extracted interactions in BioRECIPE format
	"""

	# TODO: rewrite this setting

	# rename file suffix to .xlsx
	pathname, suffix = os.path.splitext(output_file)
	if suffix != '.xlsx':
		fName = output_file.replace(suffix, ".xlsx")
	else:
		fName = output_file

	fName_stats = output_file.replace(".csv","_statements.xlsx")
	fName_VIOLIN = output_file.replace(".csv", "_statements_VIOLIN.xlsx")
	fName_FLUTE = output_file.replace(".csv", "_statements_FLUTE.xlsx")

	networkArray = np.empty((6000, 13),dtype=">U250")
	rowCount = 0
	noStatements = []

	db_refs_dict = {
		'UP': ('uniprot', 'protein'),
		'UPPRO': ('uniprot', 'protein'),
		'PF': ('pfam', 'family'),
		'CHEBI': ('chebi', 'simple-chemical'),
		'PUBCHEM': ('pubchem', 'simple-chemical'),
		'GO': ('go', 'bioprocess'),
		'MESH': ('mesh', 'bioprocess'),
		'HGNC': ('hgnc', 'gene'),
	}

	try:
		stmts = processor.statements
	except:
		raise ValueError("Processor is not supported by INDRA")

	if (stmts):
		for s in stmts:
			bs = SimpleScorer()
			j = s.to_json()
			intType = j["type"]
			bScore = bs.score_statement(s)
			valid_st = False
			try:
				# Regulator
				upstr = j["subj"]
				upstrName = upstr["name"]
				networkArray[rowCount, 0] = upstrName
				upstrDbRefs = upstr["db_refs"]

				networkArray[rowCount, 9] = bScore
				networkArray[rowCount, 8] = intType
				#networkArray[rowCount, 12] = p
				valid_st = True

				try:
					paperStats = j["evidence"]
					networkArray[rowCount, 11] = paperStats[0]["text"]

					epi_var = paperStats[0]["epistemics"]

					networkArray[rowCount, 10] = str(epi_var['direct'])

				except Exception as e:
					pass

				for db_name, db_data in db_refs_dict.items():
					try:
						upstrDb = upstrDbRefs[db_name]
						networkArray[rowCount, 1] = db_data[0]
						networkArray[rowCount, 2] = db_data[1]
						networkArray[rowCount, 3] = upstrDb
						break
					except:
						pass

				if networkArray[rowCount, 1] == 0:
					networkArray[rowCount, 3] = upstrName

			except Exception as e:
				try:
					upstr = j["sub"]
					upstrName = upstr["name"]

					networkArray[rowCount, 0] = upstrName
					upstrDbRefs = upstr["db_refs"]

					networkArray[rowCount, 9] = bScore
					networkArray[rowCount, 8] = intType
					#networkArray[rowCount, 12] = p
					valid_st = True

					try:
						paperStats = j["evidence"]
						networkArray[rowCount, 11] = paperStats[0]["text"]

						epi_var = paperStats[0]["epistemics"]

						networkArray[rowCount, 10] = str(epi_var['direct'])

					except Exception as e:
						pass

					for db_name, db_data in db_refs_dict.items():
						try:
							upstrDb = upstrDbRefs[db_name]
							networkArray[rowCount, 1] = db_data[0]
							networkArray[rowCount, 2] = db_data[1]
							networkArray[rowCount, 3] = upstrDb
							break
						except:
							pass

					if networkArray[rowCount, 1] == 0:
						networkArray[rowCount, 3] = upstrName


				except Exception as e:
					pass

			try:
				# Regulated
				downstr = j["obj"]
				downstrName = downstr["name"]
				networkArray[rowCount, 4] = downstrName
				downstrDbRefs = downstr["db_refs"]

				for db_name, db_data in db_refs_dict.items():
					try:
						downstrDb = downstrDbRefs[db_name]
						networkArray[rowCount, 5] = db_data[0]
						networkArray[rowCount, 6] = db_data[1]
						networkArray[rowCount, 7] = downstrDb
						break
					except:
						pass

				if networkArray[rowCount, 5] == 0:
					networkArray[rowCount, 7] = downstrName

			except Exception as e:
				try:
					downstr = j["enz"]
					downstrName = downstr["name"]
					networkArray[rowCount, 4] = downstrName
					downstrDbRefs = downstr["db_refs"]

					for db_name, db_data in db_refs_dict.items():
						try:
							downstrDb = downstrDbRefs[db_name]
							networkArray[rowCount, 5] = db_data[0]
							networkArray[rowCount, 6] = db_data[1]
							networkArray[rowCount, 7] = downstrDb
							break
						except:
							pass

					if networkArray[rowCount, 5] == 0:
						networkArray[rowCount, 7] = downstrName

				except Exception as e:
					pass

			if (valid_st):
				rowCount += 1

	else:
		print("No statements found")

	# h = "Regulator Name \tRegulator ID \tRegulated Name \tRegulated ID \tRegulation Type \tBelief score \tDirect? \tText"

	# original output
	fOut = pd.DataFrame(networkArray)
	fOut.columns = ["Regulator Name", "Regulator Database", "Regulator Type", "Regulator ID", "Regulated Name",
					"Regulated Database",
					"Regulated Type", "Regulated ID", "Regulation Type", "Belief score", "Connection Type", "Evidence",
					"PMCID"]

	fOut_FLUTE = pd.DataFrame(
		columns=["RegulatedName", "RegulatedDatabase", "RegulatedType", "RegulatedID", "RegulatorName",
				 "RegulatorDatabase",
				 "RegulatorID", "RegulatorType", "InteractionType", "PaperID"])

	# old VIOLIN version
	fOut_VIOLIN = pd.DataFrame(
		columns=["Element Name", "Element Type", "Database Name", "Element ID", "Location", "Location ID", "Cell Line",
				 "Cell Type", "Organism", "Positive Reg Name",
				 "Positive Reg Type", "Positive Reg ID", "Positive Reg Location", "Positive Reg Location ID",
				 "Negative Reg Name", "Negative Reg Type",
				 "Negative Reg ID", "Negative Reg Location", "Negative Reg Location ID", "Connection Type", "Mechanism",
				 "Paper ID", "Evidence"])

	# Reading output in BioRECIPE format
	fOut_BioRECIPE = pd.DataFrame(columns=BioRECIPE_reading_col)

	# convert INDRA statements to VIOLN & FLUTE formats
	for i in range(len(fOut)):
		if not fOut.loc[i, "Regulator Name"]:
			break

		fOut_FLUTE.loc[i, "RegulatedName"] = fOut.loc[i, "Regulated Name"]
		fOut_FLUTE.loc[i, "RegulatedDatabase"] = fOut.loc[i, "Regulated Database"]
		fOut_FLUTE.loc[i, "RegulatedType"] = fOut.loc[i, "Regulated Type"]
		fOut_FLUTE.loc[i, "RegulatedID"] = fOut.loc[i, "Regulated ID"]
		fOut_FLUTE.loc[i, "RegulatorName"] = fOut.loc[i, "Regulator Name"]
		fOut_FLUTE.loc[i, "RegulatorDatabase"] = fOut.loc[i, "Regulator Database"]
		fOut_FLUTE.loc[i, "RegulatorType"] = fOut.loc[i, "Regulator Type"]
		fOut_FLUTE.loc[i, "RegulatorID"] = fOut.loc[i, "Regulator ID"]

		fOut_VIOLIN.loc[i, "Element Name"] = fOut.loc[i, "Regulated Name"]
		fOut_VIOLIN.loc[i, "Element Type"] = fOut.loc[i, "Regulated Type"]
		fOut_VIOLIN.loc[i, "Element ID"] = fOut.loc[i, "Regulated ID"]

		fOut_BioRECIPE.loc[i, "Regulated Name"] = fOut.loc[i, "Regulated Name"]
		fOut_BioRECIPE.loc[i, "Regulated Type"] = fOut.loc[i, "Regulated Type"]
		fOut_BioRECIPE.loc[i, "Regulated ID"] = fOut.loc[i, "Regulated ID"]
		fOut_BioRECIPE.loc[i, "Regulated Database"] = fOut.loc[i, "Regulated Database"]
		fOut_BioRECIPE.loc[i, "Regulator Name"] = fOut.loc[i, "Regulator Name"]
		fOut_BioRECIPE.loc[i, "Regulator Type"] = fOut.loc[i, "Regulator Type"]
		fOut_BioRECIPE.loc[i, "Regulator ID"] = fOut.loc[i, "Regulator ID"]
		fOut_BioRECIPE.loc[i, "Regulator Database"] = fOut.loc[i, "Regulator Database"]

		fOut_BioRECIPE.loc[i, "Mechanism"] = fOut.loc[i, "Regulation Type"]

		if fOut.loc[i, "Regulation Type"] in ["Activation", "IncreaseAmount", "Phosphorylation"]:
			fOut_FLUTE.loc[i, "InteractionType"] = "increases"

			fOut_VIOLIN.loc[i, "Positive Reg Name"] = fOut.loc[i, "Regulator Name"]
			fOut_VIOLIN.loc[i, "Positive Reg Type"] = fOut.loc[i, "Regulator Type"]
			fOut_VIOLIN.loc[i, "Positive Reg ID"] = fOut.loc[i, "Regulator ID"]

			fOut_BioRECIPE.loc[i, "Sign"] = "positive"
		elif fOut.loc[i, "Regulation Type"] in ["Inhibition", "DecreaseAmount", "Dephosphorylation"]:
			fOut_FLUTE.loc[i, "InteractionType"] = "decreases"

			fOut_VIOLIN.loc[i, "Negative Reg Name"] = fOut.loc[i, "Regulator Name"]
			fOut_VIOLIN.loc[i, "Negative Reg Type"] = fOut.loc[i, "Regulator Type"]
			fOut_VIOLIN.loc[i, "Negative Reg ID"] = fOut.loc[i, "Regulator ID"]

			fOut_BioRECIPE.loc[i, "Sign"] = "negative"
		else:
			print("Unspecified Regulation Type: {0}".format(fOut.loc[i, "Regulation Type"]))
			fOut_FLUTE.loc[i, "InteractionType"] = "increases"

			fOut_VIOLIN.loc[i, "Positive Reg Name"] = fOut.loc[i, "Regulator Name"]
			fOut_VIOLIN.loc[i, "Positive Reg Type"] = fOut.loc[i, "Regulator Type"]
			fOut_VIOLIN.loc[i, "Positive Reg ID"] = fOut.loc[i, "Regulator ID"]

			fOut_BioRECIPE.loc[i, "Sign"] = "positive"

		if fOut.loc[i, "Connection Type"] == "True":
			fOut_VIOLIN.loc[i, "Connection Type"] = "D"
			fOut_BioRECIPE.loc[i, "Connection Type"] = "True"
		elif fOut.loc[i, "Connection Type"] == "False":
			fOut_VIOLIN.loc[i, "Connection Type"] = "I"
			fOut_BioRECIPE.loc[i, "Connection Type"] = "False"
		else:
			print("Unspecified Connection Type: {0}".format(fOut.loc[i, "Connection Type"]))
			fOut_VIOLIN.loc[i, "Connection Type"] = "I"
			fOut_BioRECIPE.loc[i, "Connection Type"] = "False"

		fOut_FLUTE.loc[i, "PaperID"] = fOut.loc[i, "PMCID"]

		fOut_VIOLIN.loc[i, "Evidence"] = fOut.loc[i, "Evidence"]
		fOut_VIOLIN.loc[i, "Paper ID"] = fOut.loc[i, "PMCID"]

		if fOut.loc[i, "Belief score"]:
			fOut_BioRECIPE.loc[i, "Score"] = fOut.loc[i, "Belief score"]
		fOut_BioRECIPE.loc[i, "Statements"] = fOut.loc[i, "Evidence"]
		fOut_BioRECIPE.loc[i, "Paper IDs"] = fOut.loc[i, "PMCID"]

	fOut_BioRECIPE.to_excel(fName, index=False)

	if indra_stats:
		fOut.to_excel(fName_stats, index=False)
	if flute:
		fOut_FLUTE = fOut_FLUTE.replace(r'^\s*$', "None", regex=True)
		fOut_FLUTE.to_excel(fName_FLUTE, index=False)
	if violin:
		fOut_VIOLIN.to_excel(fName_VIOLIN, index=False)

	print("Finished.")
	# np.savetxt(fName,networkArray,fmt="%s",encoding="utf-8",delimiter="\t", header=h,comments="")

def reach_to_biorecipe(file_name, output_file=None):

	rp = get_indra_stmts_by_reach(file_name)
	if not output_file:
		indra_stmts_reading_output(processor=rp, output_file=file_name)

def main():
	parser = argparse.ArgumentParser(description='translate INDRA statements using different systems to a reading file')
	parser.add_argument('json', type=str,
						help='A REACH json file to be read with INDRA')

	args = parser.parse_args()
	reach_to_biorecipe(args.json)

if __name__ == '__main__':
	main()
