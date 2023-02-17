import argparse
import sys
import glob
import pandas as pd
from  collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main(argv,out):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='Folder to process', default="results")
	parser.add_argument('--output', help="Output file name ['aggregated_stats.*']", default="aggregated_stats")
	parser.add_argument('--completeness_threshold', help="CheckM completeness threshold [99]", default=99)
	parser.add_argument('--contamination_threshold', help="CheckM contamination threshold [1]", default=1)
	parser.add_argument('--strain_heterogeneity_threshold', help="CheckM heterogeneity threshold [0.1]", default=0.1)
	parser.add_argument('--busco_completeness_threshold', help="BUSCO completeness threshold (%) [95]", default=95)
	parser.add_argument('--busco_duplication_threshold', help="BUSCO duplication threshold (%) [1]", default=1)
	parser.add_argument('--busco_fragmented_threshold', help="BUSCO fragmented threshold (%) [1]", default=1)
	parser.add_argument('--busco_missing_threshold', help="BUSCO missing threshold (%) [1]", default=1)
	parser.add_argument('--assembly_length_threshold_min', help="Assembly length, lower threshold (bp) [1945246]", default=1945246)
	parser.add_argument('--assembly_length_threshold_max', help="Assembly length, upper threshold (bp) [2255392]", default=2255392)
	parser.add_argument('--gc_threshold_min', help="GC content, lower threshold (%) [39.2]", default=39.2)
	parser.add_argument('--gc_threshold_max', help="GC content, upper threshold (%) [40]", default=40)
	parser.add_argument('--gap_sum_threshold', help="Total gap length (bp) [5631]", default=5631)
	parser.add_argument('--gap_count_threshold', help="Total gap count [26]", default=26)
	parser.add_argument('--perc_het_vars_threshold', help="Proportion of heterozygous variants (%) [15]", default=15)
	parser.add_argument('--scaffold_count_threshold', help="Maximum number of scaffolds per assembly [286]", default=286)
	parser.add_argument('--scaffold_N50_threshold', help="Minimum scaffold N50 [24454]", default=24454)
	parser.add_argument('--MASH_hit', help="Closest MASH hit (top 5) [Streptococcus pneumoniae]", default="Streptococcus pneumoniae")
	args = parser.parse_args(argv)
	result = add_label_column(args.input, args.completeness_threshold, args.contamination_threshold, args.strain_heterogeneity_threshold, args.busco_completeness_threshold, args.busco_duplication_threshold, args.busco_fragmented_threshold, args.busco_missing_threshold, args.assembly_length_threshold_min, args.assembly_length_threshold_max, args.gc_threshold_min, args.gc_threshold_max, args.gap_sum_threshold, args.gap_count_threshold, args.perc_het_vars_threshold, args.scaffold_count_threshold, args.scaffold_N50_threshold, args.MASH_hit)
	dfx = pd.DataFrame([sub.split(",") for sub in result])
	dfx.columns =['sample', 'metric', 'result', 'status']
	dfx.to_csv(args.output+".long.txt", sep=',', index=False)
	rx3 = aggregate(dfx)
	rx3.to_csv(args.output+".wide.txt", sep=',', index=False)
	plot_all(rx3)

def add_label_column(input, completeness_threshold, contamination_threshold, strain_heterogeneity_threshold, busco_completeness_threshold, busco_duplication_threshold, busco_fragmented_threshold, busco_missing_threshold, assembly_length_threshold_min, assembly_length_threshold_max, gc_threshold_min, gc_threshold_max, gap_sum_threshold, gap_count_threshold, perc_het_vars_threshold, scaffold_count_threshold, scaffold_N50_threshold, MASH_hit):
	input2 = input
	pattern = f"{input2}/*/*_results.txt"
	file_paths = glob.glob(pattern)
	contents = ""
	for path in file_paths:
		with open(path, "r") as f:
			contents += f.read()
	lines = contents.splitlines()
	result = []
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
		elif metric == "CHECKM_Strain_heterogeneity":
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
			if str(value) == MASH_hit:
				result.append(f"{sample_id},{metric},{value},PASS")
			else:
				result.append(f"{sample_id},{metric},{value},FAIL")
		else:
			result.append(f"{sample_id},{metric},{value},NA")
	return result

def aggregate(input):
	df = input
	df2 = df.pivot(index='sample', columns='metric', values='result')
	xf1 = df.loc[df['status']=='FAIL'].groupby('sample')['status'].agg({'count'})
	df3 = df2.merge(xf1, left_on='sample', right_on='sample',how='left').fillna("NA")
	df3[['count']] = df3[['count']].replace("NA", 0)
	df3 = df3.rename(columns={'count': 'fail_counts'})
	return df3

def plot_all(input):
	df=input
	sns.set_style("white")
	sns.set_theme()
	sns.set(style="ticks")

	pdf_pages = PdfPages("aggregated_plots.pdf")
	fig, ((axs1, axs2), (axs3, axs4), (axs5, axs6), (axs7, axs8)  ) = plt.subplots(4, 2, figsize=(8, 12))

	if "assembly_length_bp" in df.columns:
		sns.histplot(data=df, x="assembly_length_bp", binwidth=10000, ax=axs1)
		axs1.set_xlabel("Assembly length (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs1.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs1.axvline(x=1945246, color="grey", lw=1, linestyle="dashed")
		axs1.axvline(x=2255392, color="grey", lw=1, linestyle="dashed")
		axs1.tick_params(axis='x', labelsize=6)
		axs1.tick_params(axis='y', labelsize=6)
	else:
		axs1.set_visible(False)

	if "GC_perc" in df.columns:
		sns.histplot(data=df, x="GC_perc", binwidth=0.1, ax=axs2)
		axs2.set_xlabel("GC content (%)", fontsize=8, fontdict={"weight": "bold"})
		axs2.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs2.axvline(x=40, color="grey", lw=1, linestyle="dashed")
		axs2.axvline(x=39.2, color="grey", lw=1, linestyle="dashed")
		axs2.tick_params(axis='x', labelsize=6)
		axs2.tick_params(axis='y', labelsize=6)
	else:
		axs2.set_visible(False)

	if "contig_count" in df.columns:
		sns.histplot(data=df, x="contig_count", binwidth=10, ax=axs3)
		axs3.set_xlabel("Contigs (count)", fontsize=8, fontdict={"weight": "bold"})
		axs3.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs3.tick_params(axis='x', labelsize=6)
		axs3.tick_params(axis='y', labelsize=6)
	else:
		axs3.set_visible(False)

	if "scaffold_count" in df.columns:
		sns.histplot(data=df, x="scaffold_count", binwidth=10, ax=axs4)
		axs4.set_xlabel("Scaffold (count)", fontsize=8, fontdict={"weight": "bold"})
		axs4.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs4.axvline(x=286, color="grey", lw=1, linestyle="dashed")
		axs4.tick_params(axis='x', labelsize=6)
		axs4.tick_params(axis='y', labelsize=6)
	else:
		axs4.set_visible(False)

	if "contig_N50_bp" in df.columns:
		sns.histplot(data=df, x="contig_N50_bp", binwidth=10000, ax=axs5)
		axs5.set_xlabel("Contig N50 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs5.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs5.tick_params(axis='x', labelsize=6)
		axs5.tick_params(axis='y', labelsize=6)
	else:
		axs5.set_visible(False)

	if "scaffold_N50_bp" in df.columns:
		sns.histplot(data=df, x="scaffold_N50_bp", binwidth=10000, ax=axs6)
		axs6.set_xlabel("Scaffold N50 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs6.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs6.axvline(x=24454, color="grey", lw=1, linestyle="dashed")
		axs6.tick_params(axis='x', labelsize=6)
		axs6.tick_params(axis='y', labelsize=6)
	else:
		axs6.set_visible(False)

	if "contig_N90_bp" in df.columns:
		sns.histplot(data=df, x="contig_N90_bp", binwidth=10000, ax=axs7)
		axs7.set_xlabel("Contig N90 (bp)", fontsize=8, fontdict={"weight": "bold"})
		axs7.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		axs7.tick_params(axis='x', labelsize=6)
		axs7.tick_params(axis='y', labelsize=6)
	else:
		axs7.set_visible(False)

	if "scaffold_N90_bp" in df.columns:
		sns.histplot(data=df, x="scaffold_N90_bp", binwidth=10000, ax=axs8)
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

	fig, ((bsx7,bsx8), (bsx1, bsx2), (bsx3, bsx4), (bsx5, bsx6)  ) = plt.subplots(4, 2, figsize=(8, 12))

	if "gaps_count" in df.columns:
		sns.histplot(data=df, x="gaps_count", binwidth=0.1, ax=bsx7)
		bsx7.set_xlabel("Gaps (count)", fontsize=8, fontdict={"weight": "bold"})
		bsx7.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx7.axvline(x=26, color="grey", lw=1, linestyle="dashed")
		bsx7.tick_params(axis='x', labelsize=6)
		bsx7.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx7.set_visible(False)

	if "gaps_sum_bp" in df.columns:
		sns.histplot(data=df, x="gaps_sum_bp", binwidth=1000, ax=bsx8)
		bsx8.set_xlabel("Gap length (bp)", fontsize=8, fontdict={"weight": "bold"})
		bsx8.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx8.axvline(x=5631, color="grey", lw=1, linestyle="dashed")
		bsx8.tick_params(axis='x', labelsize=6)
		bsx8.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx8.set_visible(False)

	if "BUSCO_complete_single_copy" in df.columns:
		sns.histplot(data=df, x="BUSCO_complete_single_copy", binwidth=0.1, ax=bsx1)
		bsx1.set_xlabel("BUSCO: Complete and single-copy (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx1.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx1.axvline(x=95, color="grey", lw=1, linestyle="dashed")
		bsx1.tick_params(axis='x', labelsize=6)
		bsx1.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx1.set_visible(False)

	if "BUSCO_complete_duplicated" in df.columns:
		sns.histplot(data=df, x="BUSCO_complete_duplicated", binwidth=0.1, ax=bsx2)
		bsx2.set_xlabel("BUSCO: Complete and duplicated (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx2.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx2.axvline(x=1, color="grey", lw=1, linestyle="dashed")
		bsx2.tick_params(axis='x', labelsize=6)
		bsx2.tick_params(axis='y', labelsize=6)
	else:
		bsx2.set_visible(False)

	if "BUSCO_fragmented" in df.columns:
		sns.histplot(data=df, x="BUSCO_fragmented", binwidth=0.1, ax=bsx3)
		bsx3.set_xlabel("BUSCO: Fragmented (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx3.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx3.axvline(x=1, color="grey", lw=1, linestyle="dashed")
		bsx3.tick_params(axis='x', labelsize=6)
		bsx3.tick_params(axis='y', labelsize=6)
		sns.despine()
	else:
		bsx3.set_visible(False)

	if "BUSCO_missing" in df.columns:
		sns.histplot(data=df, x="BUSCO_missing", binwidth=0.1, ax=bsx4)
		bsx4.set_xlabel("BUSCO: Missing (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx4.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx4.axvline(x=1, color="grey", lw=1, linestyle="dashed")
		bsx4.tick_params(axis='x', labelsize=6)
		bsx4.tick_params(axis='y', labelsize=6)
	else:
		bsx4.set_visible(False)

	if "perc_het_vars" in df.columns:
		sns.histplot(data=df, x="perc_het_vars", binwidth=0.1, ax=bsx5)
		bsx5.set_xlabel("Proportion of heterozygous variants (%)", fontsize=8, fontdict={"weight": "bold"})
		bsx5.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx5.axvline(x=10, color="grey", lw=1, linestyle="dashed")
		bsx5.tick_params(axis='x', labelsize=6)
		bsx5.tick_params(axis='y', labelsize=6)
	else:
		bsx5.set_visible(False)

	if "MASH_hit" in df.columns:
		counts = df["MASH_hit"].value_counts().reset_index()
		counts.columns = ['MASH_hit', 'countz']
		sns.barplot(counts, x="MASH_hit", y='countz', ax=bsx6)
		bsx6.set_xlabel("Top 5 MASH hits", fontsize=8, fontdict={"weight": "bold"})
		bsx6.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		bsx6.tick_params(axis='x', labelsize=6, rotation=45)
		labels = bsx6.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		bsx6.set_xticklabels(labels)
		bsx6.tick_params(axis='y', labelsize=6)
	else:
		bsx6.set_visible(False)

	sns.despine()
	fig.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9, wspace=0.1)
	fig.subplots_adjust(hspace=0.4, wspace =0.4)
	pdf_pages.savefig(fig)

	fig, ((csv1, csv2),(csv5,_), (csv3,csv4)) = plt.subplots(3, 2, figsize=(8, 9))

	if "pneumoKITy_serotype" in df.columns:
		top_serotypes = df["pneumoKITy_serotype"].value_counts().nlargest(10).index.tolist()
		df_top = df[df["pneumoKITy_serotype"].isin(top_serotypes)]
		counts = df_top["pneumoKITy_serotype"].value_counts().reset_index()
		counts.columns = ['pneumoKITy_serotype', 'countz']
		sns.barplot(counts, x="pneumoKITy_serotype", y='countz', ax=csv1)
		csv1.set_xlabel("PnuemoKITy serotype", fontsize=8, fontdict={"weight": "bold"})
		csv1.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv1.tick_params(axis='x', labelsize=6, rotation=45)
		labels = csv1.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		csv1.set_xticklabels(labels)
		csv1.tick_params(axis='y', labelsize=6)
	else:
		csv1.set_visible(False)

	if "seroBA_serotype" in df.columns:
		top_serotypes = df["seroBA_serotype"].value_counts().nlargest(10).index.tolist()
		df_top = df[df["seroBA_serotype"].isin(top_serotypes)]
		counts = df_top["seroBA_serotype"].value_counts().reset_index()
		counts.columns = ['seroBA_serotype', 'countz']
		sns.barplot(counts, x="seroBA_serotype", y='countz', ax=csv2)
		csv2.set_xlabel("seroBA serotype", fontsize=8, fontdict={"weight": "bold"})
		csv2.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv2.tick_params(axis='x', labelsize=6, rotation=45)
		labels = csv2.get_xticklabels()
		for label in labels:
			label.set_ha('right')
		csv2.set_xticklabels(labels)
		csv2.tick_params(axis='y', labelsize=6)
	else:
		csv2.set_visible(False)

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
		sns.histplot(data=df, x="n_genes", binwidth=1, ax=csv3)
		csv3.set_xlabel("Predicted genes (count)", fontsize=8, fontdict={"weight": "bold"})
		csv3.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv3.tick_params(axis='x', labelsize=6)
		csv3.tick_params(axis='y', labelsize=6)
	else:
		csv3.set_visible(False)

	if "n_genes" in df.columns:
		sns.histplot(data=df, x="fail_counts", binwidth=1, ax=csv4)
		csv4.set_xlabel("Failures per sample", fontsize=8, fontdict={"weight": "bold"})
		csv4.set_ylabel("Frequency", fontsize=8, fontdict={"weight": "bold"})
		csv4.tick_params(axis='x', labelsize=6)
		csv4.tick_params(axis='y', labelsize=6)
	else:
		csv4.set_visible(False)

	_.set_visible(False)
	_.spines["left"].set_visible(False)
	_.spines["right"].set_visible(False)
	_.spines["top"].set_visible(False)
	_.spines["bottom"].set_visible(False)
	_.xaxis.set_visible(False)
	_.yaxis.set_visible(False)
	
	sns.despine()
	fig.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9, wspace=0.1)
	fig.subplots_adjust(hspace=0.4, wspace =0.4)
	pdf_pages.savefig(fig)
	pdf_pages.close()

if __name__ == '__main__':
	main(sys.argv[1:], sys.stdout)
