import argparse
import sys
import glob
import subprocess
import pandas as pd
from  collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main(argv,out): 
	# Add input arguments
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='Folder to process', default="results")
	parser.add_argument('--output', help="Output file name ['aggregated_stats.*']", default="aggregated_stats")
	parser.add_argument('--completeness_threshold', help="CheckM completeness threshold [90]", default=90)
	parser.add_argument('--contamination_threshold', help="CheckM contamination threshold [5]", default=5)
	parser.add_argument('--strain_heterogeneity_threshold', help="CheckM heterogeneity threshold [40]", default=40)
	parser.add_argument('--busco_completeness_threshold', help="BUSCO completeness threshold (%) [95]", default=95)
	parser.add_argument('--busco_duplication_threshold', help="BUSCO duplication threshold (%) [1]", default=1)
	parser.add_argument('--busco_fragmented_threshold', help="BUSCO fragmented threshold (%) [1]", default=1)
	parser.add_argument('--busco_missing_threshold', help="BUSCO missing threshold (%) [1]", default=1)
	parser.add_argument('--assembly_length_threshold_min', help="Assembly length, lower threshold (bp) [1945246]", default=1945246)
	parser.add_argument('--assembly_length_threshold_max', help="Assembly length, upper threshold (bp) [2255392]", default=2255392)
	parser.add_argument('--gc_threshold_min', help="GC content, lower threshold (%) [39.2]", default=39.2)
	parser.add_argument('--gc_threshold_max', help="GC content, upper threshold (%) [40]", default=40)
	parser.add_argument('--gap_sum_threshold', help="Total gap length (bp) [5631]", default=5631)
	parser.add_argument('--gap_count_threshold', help="Total gap count [50]", default=50)
	parser.add_argument('--perc_het_vars_threshold', help="Proportion of heterozygous variants (%) [15]", default=15)
	parser.add_argument('--scaffold_count_threshold', help="Maximum number of scaffolds per assembly [286]", default=286)
	parser.add_argument('--scaffold_N50_threshold', help="Minimum scaffold N50 [24454]", default=24454)
	parser.add_argument('--target_species', help="Closest MASH hit (top 5) [Streptococcus pneumoniae]", default="Streptococcus pneumoniae")
	parser.add_argument('--target_genus', help="Closest MASH hit (top 5) [Streptococcus]", default="Streptococcus")
	parser.add_argument('--kraken_match_threshold_species', help="Minimum percentage of reads coverged by the --target_species [25]", default=25)
	parser.add_argument('--kraken_match_threshold_genus', help="Minimum percentage of reads coverged by the --target_genus [50]", default=50)

	args = parser.parse_args(argv)
	# Make a directory for the aggregated output files
	instring="mkdir -p "+args.input+"/"+args.output+"/"
	subprocess.call(instring, shell=True)
	# Write a command and run multiqc
	mqc = "multiqc --outdir "+args.input+"/"+args.output+"/ "+args.input+"/*/"
	subprocess.call(mqc, shell=True)
	# Run the first function, to parse output files, aggregate results, identify failed assemblies
	result = add_label_column(args.input, args.completeness_threshold, args.contamination_threshold, args.strain_heterogeneity_threshold, args.busco_completeness_threshold, args.busco_duplication_threshold, args.busco_fragmented_threshold, args.busco_missing_threshold, args.assembly_length_threshold_min, args.assembly_length_threshold_max, args.gc_threshold_min, args.gc_threshold_max, args.gap_sum_threshold, args.gap_count_threshold, args.perc_het_vars_threshold, args.scaffold_count_threshold, args.scaffold_N50_threshold, args.target_species, args.kraken_match_threshold_species, args.kraken_match_threshold_genus, args.target_genus)
	# Reformat output
	dfx = pd.DataFrame([sub.split(",") for sub in result])
	dfx.columns =['sample', 'metric', 'result', 'status']
	# Write to file
	dfx.to_csv(args.input+"/"+args.output+"/"+args.output+".long.csv", sep=',', index=False)
	# Aggregate results to a reformatted table and write to file
	rx3 = aggregate(dfx)
	rx3.to_csv(args.input+"/"+args.output+"/"+args.output+".wide.csv", sep=',', index=True)
	# Make plots
	plot_all(args.input+"/"+args.output+"/"+args.output+".pdf", args.input+"/"+args.output+"/"+args.output+".wide.csv",rx3, args.completeness_threshold, args.contamination_threshold, args.strain_heterogeneity_threshold, args.busco_completeness_threshold, args.busco_duplication_threshold, args.busco_fragmented_threshold, args.busco_missing_threshold, args.assembly_length_threshold_min, args.assembly_length_threshold_max, args.gc_threshold_min, args.gc_threshold_max, args.gap_sum_threshold, args.gap_count_threshold, args.perc_het_vars_threshold, args.scaffold_count_threshold, args.scaffold_N50_threshold)

def add_label_column(input, completeness_threshold, contamination_threshold, strain_heterogeneity_threshold, busco_completeness_threshold, busco_duplication_threshold, busco_fragmented_threshold, busco_missing_threshold, assembly_length_threshold_min, assembly_length_threshold_max, gc_threshold_min, gc_threshold_max, gap_sum_threshold, gap_count_threshold, perc_het_vars_threshold, scaffold_count_threshold, scaffold_N50_threshold, target_species, kraken_match_threshold_species, kraken_match_threshold_genus, target_genus):
	# Organise the importatation of data from multiple folders
	pattern = f"{input}/*/*_results.txt"
	file_paths = glob.glob(pattern)
	contents = ""
	# Read data in and split lines
	for path in file_paths:
		with open(path, "r") as f:
			contents += f.read()
	lines = contents.splitlines()
	result = []
	# For each relevant metric, identify whether it's a PASS or FAIL based on default or user defined thresholds
	for line in lines:
		columns = line.strip().split(',')
		sample_id, metric, value, backup = columns[0], columns[1], columns[2], columns[3]
		if metric == "CHECKM_Completeness":
			if float(value) >= completeness_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "CHECKM_Contamination":
			if float(value) <= contamination_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "v":
			if float(value) <= strain_heterogeneity_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "BUSCO_complete_single_copy":
			if float(value) >= busco_completeness_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "BUSCO_complete_duplicated":
			if float(value) <= busco_duplication_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "BUSCO_fragmented":
			if float(value) <= busco_fragmented_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "BUSCO_missing":
			if float(value) <= busco_missing_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "pneumoKITy_serotype":
			if str(backup) == "PASS":
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "assembly_length_bp":
			if assembly_length_threshold_min  <= float(value) <= assembly_length_threshold_max :
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "GC_perc":
			if gc_threshold_min  <= float(value) <= gc_threshold_max :
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "gaps_sum_bp":
			if float(value) <= gap_sum_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "gaps_count":
			if float(value) <= gap_count_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "perc_het_vars":
			if float(value) <= perc_het_vars_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "scaffold_count":
			if float(value) <= scaffold_count_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "scaffold_N50_bp":
			if float(value) >= scaffold_N50_threshold:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "MASH_hit":
			if str(value) == target_species:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		elif metric == "kraken2_species":
			colx = value.strip().split(':')
			species,perc = colx[0], colx[1]
			if str(species) == target_species and float(perc)>=kraken_match_threshold_species :
				result.append(f"{sample_id},{metric},{species},PASS")
			else:
				result.append(f"{sample_id},{metric},{species},FAIL")
		elif metric == "kraken2_genus":
			colx = value.strip().split(':')
			species,perc = colx[0], colx[1]
			if str(species) == target_genus and float(perc)>=kraken_match_threshold_genus :
				result.append(f"{sample_id},{metric},{species},PASS")
			else:
				result.append(f"{sample_id},{metric},{species},FAIL")
		elif "mlst_allele" in metric:
			if ';' in str(value) or '?' in str(value) or str(value).strip() == '':
				result.append(f"{sample_id},{metric},{value},FAIL")
			else:
				result.append(f"{sample_id},{metric},{value},PASS")
		else:
			result.append(f"{sample_id},{metric},{value},NA")
	return result

def aggregate(input):
	# Merge data table from add_label_column with a table of counts of FAIL metrics (per-sample)
	df = input
	df2 = df.pivot(index='sample', columns='metric', values='result')
	xf1 = df.loc[df['status']=='FAIL'].groupby('sample')['status'].agg({'count'})
	df3 = df2.merge(xf1, left_on='sample', right_on='sample',how='left').fillna("NA")
	df3[['count']] = df3[['count']].replace("NA", 0)
	df3 = df3.rename(columns={'count': 'fail_counts'})
	return df3

def plot_all(outputstring, inputstring,input, completeness_threshold, contamination_threshold, strain_heterogeneity_threshold, busco_completeness_threshold, busco_duplication_threshold, busco_fragmented_threshold, busco_missing_threshold, assembly_length_threshold_min, assembly_length_threshold_max, gc_threshold_min, gc_threshold_max, gap_sum_threshold, gap_count_threshold, perc_het_vars_threshold, scaffold_count_threshold, scaffold_N50_threshold):
	df=pd.read_csv(inputstring)
	sns.set_style("white")
	sns.set_theme()
	sns.set(style="ticks")
	# Open PDF and plot figures across 3 pages (up to 8 plots per page)
	pdf_pages = PdfPages(outputstring)
	fig, ((axs1, axs2), (axs3, axs4), (axs5, axs6), (axs7, axs8)  ) = plt.subplots(4, 2, figsize=(8, 12))

	if "assembly_length_bp" in df.columns:
		sns.histplot(data=df, x="assembly_length_bp", bins=50, ax=axs1)
		axs1.set_xlabel("Assembly length (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs1.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs1.axvline(x=assembly_length_threshold_min, color="grey", lw=1, linestyle="dashed")
		axs1.axvline(x=assembly_length_threshold_max, color="grey", lw=1, linestyle="dashed")
		axs1.tick_params(axis='x', labelsize=6)
		axs1.tick_params(axis='y', labelsize=6)
	else:
		axs1.set_visible(False)

	if "GC_perc" in df.columns:
		sns.histplot(data=df, x="GC_perc", bins=50, ax=axs2)
		axs2.set_xlabel("GC content (%)", fontsize=8, fontdict={"weight": "bold"})
		axs2.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs2.axvline(x=gc_threshold_max, color="grey", lw=1, linestyle="dashed")
		axs2.axvline(x=gc_threshold_min, color="grey", lw=1, linestyle="dashed")
		axs2.tick_params(axis='x', labelsize=6)
		axs2.tick_params(axis='y', labelsize=6)
	else:
		axs2.set_visible(False)

	if "contig_count" in df.columns:
		sns.histplot(data=df, x="contig_count", bins=50, ax=axs3)
		axs3.set_xlabel("Contigs (count)", fontsize=8, fontdict={"weight": "bold"})
		axs3.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs3.tick_params(axis='x', labelsize=6)
		axs3.tick_params(axis='y', labelsize=6)
	else:
		axs3.set_visible(False)

	if "scaffold_count" in df.columns:
		sns.histplot(data=df, x="scaffold_count", bins=50, ax=axs4)
		axs4.set_xlabel("Scaffold (count)", fontsize=8, fontdict={"weight": "bold"})
		axs4.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs4.axvline(x=scaffold_count_threshold, color="grey", lw=1, linestyle="dashed")
		axs4.tick_params(axis='x', labelsize=6)
		axs4.tick_params(axis='y', labelsize=6)
	else:
		axs4.set_visible(False)

	if "contig_N50_bp" in df.columns:
		sns.histplot(data=df, x="contig_N50_bp", bins=50, ax=axs5)
		axs5.set_xlabel("Contig N50 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs5.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs5.tick_params(axis='x', labelsize=6)
		axs5.tick_params(axis='y', labelsize=6)
	else:
		axs5.set_visible(False)

	if "scaffold_N50_bp" in df.columns:
		sns.histplot(data=df, x="scaffold_N50_bp", bins=50, ax=axs6)
		axs6.set_xlabel("Scaffold N50 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs6.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs6.axvline(x=scaffold_N50_threshold, color="grey", lw=1, linestyle="dashed")
		axs6.tick_params(axis='x', labelsize=6)
		axs6.tick_params(axis='y', labelsize=6)
	else:
		axs6.set_visible(False)

	if "contig_N90_bp" in df.columns:
		sns.histplot(data=df, x="contig_N90_bp", bins=50, ax=axs7)
		axs7.set_xlabel("Contig N90 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs7.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs7.tick_params(axis='x', labelsize=6)
		axs7.tick_params(axis='y', labelsize=6)
	else:
		axs7.set_visible(False)

	if "scaffold_N90_bp" in df.columns:
		sns.histplot(data=df, x="scaffold_N90_bp", bins=50, ax=axs8)
		axs8.set_xlabel("Scaffold N90 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs8.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs8.tick_params(axis='x', labelsize=6)
		axs8.tick_params(axis='y', labelsize=6)
	else:
		axs8.set_visible(False)

	sns.despine()
	fig.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9, wspace=0.1)
	fig.subplots_adjust(hspace=0.4, wspace =0.4)
	pdf_pages.savefig(fig)

	fig, ((bsx7,bsx8), (bsx1, bsx2), (bsx3, bsx4), (bsx9, bsx10)  ) = plt.subplots(4, 2, figsize=(8, 12))

	if "gaps_count" in df.columns:
		sns.histplot(data=df, x="gaps_count", bins=50, ax=bsx7)
		bsx7.set_xlabel("Gaps (count)", fontsize=8, fontdict={"weight": "bold"})
		bsx7.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx7.axvline(x=gap_count_threshold, color="grey", lw=1, linestyle="dashed")
		bsx7.tick_params(axis='x', labelsize=6)
		bsx7.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx7.set_visible(False)

	if "gaps_sum_bp" in df.columns:
		sns.histplot(data=df, x="gaps_sum_bp", bins=50, ax=bsx8)
		bsx8.set_xlabel("Gap length (bp)", fontsize=8, fontdict={"weight": "bold"})
		bsx8.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx8.axvline(x=gap_sum_threshold, color="grey", lw=1, linestyle="dashed")
		bsx8.tick_params(axis='x', labelsize=6)
		bsx8.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx8.set_visible(False)

	if "BUSCO_complete_single_copy" in df.columns:
		sns.histplot(data=df, x="BUSCO_complete_single_copy", bins=30, ax=bsx1)
		bsx1.set_xlabel("BUSCO: Complete and single-copy (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx1.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx1.axvline(x=busco_completeness_threshold, color="grey", lw=1, linestyle="dashed")
		bsx1.tick_params(axis='x', labelsize=6)
		bsx1.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx1.set_visible(False)

	if "BUSCO_complete_duplicated" in df.columns:
		sns.histplot(data=df, x="BUSCO_complete_duplicated", bins=30, ax=bsx2)
		bsx2.set_xlabel("BUSCO: Complete and duplicated (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx2.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx2.axvline(x=busco_duplication_threshold, color="grey", lw=1, linestyle="dashed")
		bsx2.tick_params(axis='x', labelsize=6)
		bsx2.tick_params(axis='y', labelsize=6)
	else:
		bsx2.set_visible(False)

	if "BUSCO_fragmented" in df.columns:
		sns.histplot(data=df, x="BUSCO_fragmented", bins=30, ax=bsx3)
		bsx3.set_xlabel("BUSCO: Fragmented (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx3.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx3.axvline(x=busco_fragmented_threshold, color="grey", lw=1, linestyle="dashed")
		bsx3.tick_params(axis='x', labelsize=6)
		bsx3.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx3.set_visible(False)

	if "BUSCO_missing" in df.columns:
		sns.histplot(data=df, x="BUSCO_missing", bins=30, ax=bsx4)
		bsx4.set_xlabel("BUSCO: Missing (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx4.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx4.axvline(x=busco_missing_threshold, color="grey", lw=1, linestyle="dashed")
		bsx4.tick_params(axis='x', labelsize=6)
		bsx4.tick_params(axis='y', labelsize=6)
	else:
		bsx4.set_visible(False)

	if "CHECKM_Completeness" in df.columns:
		sns.histplot(data=df, x="CHECKM_Completeness", bins=30, ax=bsx10)
		bsx10.set_xlabel("CHECKM completeness (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx10.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx10.axvline(x=completeness_threshold, color="grey", lw=1, linestyle="dashed")
		bsx10.tick_params(axis='x', labelsize=6)
		bsx10.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx10.set_visible(False)

	if "CHECKM_Contamination" in df.columns:
		sns.histplot(data=df, x="CHECKM_Contamination", bins=30, ax=bsx9)
		bsx9.set_xlabel("CHECKM contamination (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx9.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx9.axvline(x=contamination_threshold, color="grey", lw=1, linestyle="dashed")
		bsx9.tick_params(axis='x', labelsize=6)
		bsx9.tick_params(axis='y', labelsize=6)
	else:
		bsx9.set_visible(False)

	# Remove extra axes
	sns.despine()
	# Fix spacing on page and between plots
	fig.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9, wspace=0.1)
	fig.subplots_adjust(hspace=0.4, wspace =0.4)
	pdf_pages.savefig(fig)

	fig, ((csx12,csx7), (csv3,csv4),(csv5,csx14),(csx13,csx8)) = plt.subplots(4, 2, figsize=(8, 12))

	if "CHECKM_Strain_heterogeneity" in df.columns:
		sns.histplot(data=df, x="CHECKM_Strain_heterogeneity", bins=30, ax=csx12)
		csx12.set_xlabel("CHECKM strain heterogeneity (%)", fontsize=8, fontdict={"weight": "bold"})
		csx12.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csx12.axvline(x=strain_heterogeneity_threshold, color="grey", lw=1, linestyle="dashed")
		csx12.tick_params(axis='x', labelsize=6)
		csx12.tick_params(axis='y', labelsize=6)
	else:
		csx12.set_visible(False)

	if "GPSC" in df.columns:
		top_GPSC = df["GPSC"].value_counts().nlargest(10).index.tolist()
		df_top = df[df["GPSC"].isin(top_GPSC)]
		counts = df_top["GPSC"].value_counts().reset_index()
		counts.columns = ['GPSC', 'countz']
		sns.barplot(counts, x="GPSC", y='countz', ax=csv5)
		csv5.set_xlabel("GPSC", fontsize=8, fontdict={"weight": "bold"})
		csv5.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv5.tick_params(axis='x', labelsize=6, rotation=45)
		labels = csv5.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		csv5.set_xticklabels(labels)
		csv5.tick_params(axis='y', labelsize=6)
	else:
		csv5.set_visible(False)

	if "n_genes" in df.columns:
		sns.histplot(data=df, x="n_genes", bins=50, ax=csv3)
		csv3.set_xlabel("Predicted genes (count)", fontsize=8, fontdict={"weight": "bold"})
		csv3.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv3.tick_params(axis='x', labelsize=6)
		csv3.tick_params(axis='y', labelsize=6)
	else:
		csv3.set_visible(False)

	if "fail_counts" in df.columns:
		sns.histplot(data=df, x="fail_counts", bins=50, ax=csv4)
		csv4.set_xlabel("Failures per sample", fontsize=8, fontdict={"weight": "bold"})
		csv4.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv4.tick_params(axis='x', labelsize=6)
		csv4.tick_params(axis='y', labelsize=6)
	else:
		csv4.set_visible(False)

	if "perc_het_vars" in df.columns:
		sns.histplot(data=df, x="perc_het_vars", bins=30, ax=csx7)
		csx7.set_xlabel("Proportion of heterozygous variants (%)", fontsize=8, fontdict={"weight": "bold"})
		csx7.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csx7.axvline(x=perc_het_vars_threshold, color="grey", lw=1, linestyle="dashed")
		csx7.tick_params(axis='x', labelsize=6)
		csx7.tick_params(axis='y', labelsize=6)
	else:
		csx7.set_visible(False)

	if "MASH_hit" in df.columns:
		counts = df["MASH_hit"].value_counts().reset_index()
		counts.columns = ['MASH_hit', 'countz']
		sns.barplot(counts, x="MASH_hit", y='countz', ax=csx8)
		csx8.set_xlabel("Top 5 MASH hits", fontsize=8, fontdict={"weight": "bold"})
		csx8.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csx8.tick_params(axis='x', labelsize=5, rotation=45)
		labels = csx8.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		csx8.set_xticklabels(labels)
		csx8.tick_params(axis='y', labelsize=6)
	else:
		csx8.set_visible(False)

	if "kraken2_species" in df.columns:
		counts = df["kraken2_species"].value_counts().reset_index()
		counts.columns = ['kraken2_species', 'countz']
		sns.barplot(counts, x="kraken2_species", y='countz', ax=csx13)
		csx13.set_xlabel("Top Kraken2 species match", fontsize=8, fontdict={"weight": "bold"})
		csx13.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csx13.tick_params(axis='x', labelsize=5, rotation=45)
		labels = csx13.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		csx13.set_xticklabels(labels)
		csx13.tick_params(axis='y', labelsize=6)
	else:
		csx13.set_visible(False)

	if "kraken2_genus" in df.columns:
		counts = df["kraken2_genus"].value_counts().reset_index()
		counts.columns = ['kraken2_genus', 'countz']
		sns.barplot(counts, x="kraken2_genus", y='countz', ax=csx14)
		csx14.set_xlabel("Top Kraken2 genus match", fontsize=8, fontdict={"weight": "bold"})
		csx14.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csx14.tick_params(axis='x', labelsize=5, rotation=45)
		labels = csx14.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		csx14.set_xticklabels(labels)
		csx14.tick_params(axis='y', labelsize=6)
	else:
		csx14.set_visible(False)
	# Remove spacing
	sns.despine()
	# Fix spacing on page and between plots
	fig.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9, wspace=0.1)
	fig.subplots_adjust(hspace=0.8, wspace =0.4)
	pdf_pages.savefig(fig)

	fig, (dsx1,dsx2) = plt.subplots(1, 2, figsize=(8, 3))


	if "pneumoKITy_serotype" in df.columns:
		top_serotypes = df["pneumoKITy_serotype"].value_counts().nlargest(10).index.tolist()
		df_top = df[df["pneumoKITy_serotype"].isin(top_serotypes)]
		counts = df_top["pneumoKITy_serotype"].value_counts().reset_index()
		counts.columns = ['pneumoKITy_serotype', 'countz']
		sns.barplot(counts, x="pneumoKITy_serotype", y='countz', ax=dsx1)
		dsx1.set_xlabel("PneumoKITy serotype", fontsize=8, fontdict={"weight": "bold"})
		dsx1.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		dsx1.tick_params(axis='x', labelsize=6, rotation=45)
		labels = dsx1.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		dsx1.set_xticklabels(labels)
		dsx1.tick_params(axis='y', labelsize=6)
	else:
		dsx1.set_visible(False)

	if "seroBA_serotype" in df.columns:
		top_serotypes = df["seroBA_serotype"].value_counts().nlargest(10).index.tolist()
		df_top = df[df["seroBA_serotype"].isin(top_serotypes)]
		counts = df_top["seroBA_serotype"].value_counts().reset_index()
		counts.columns = ['seroBA_serotype', 'countz']
		sns.barplot(counts, x="seroBA_serotype", y='countz', ax=dsx2)
		dsx2.set_xlabel("seroBA serotype", fontsize=8, fontdict={"weight": "bold"})
		dsx2.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		dsx2.tick_params(axis='x', labelsize=6, rotation=45)
		labels = dsx2.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		dsx2.set_xticklabels(labels)
		dsx2.tick_params(axis='y', labelsize=6)
	else:
		dsx2.set_visible(False)

	sns.despine()
	# Fix spacing on page and between plots
	fig.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9, wspace=0.1)
	fig.subplots_adjust(hspace=0.8, wspace =0.4)
	pdf_pages.savefig(fig)

	# Close PDF
	pdf_pages.close()

if __name__ == '__main__':
	main(sys.argv[1:], sys.stdout)
