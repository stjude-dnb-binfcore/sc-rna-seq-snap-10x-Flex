#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author:
	Name: 			Antonia Chroni
	Email: 			antonia.chroni@stjude.org
	Affiliation: 	St. Jude Children's Research Hospital, Memphis, TN
	Date: 			May 27, 2025 
"""

import os, sys, argparse, glob, numpy, pandas
import pandas as pd


def dir_path(string):
	if os.path.isdir(string):
		return string
	else:
		raise NotADirectoryError(string)


# Interpreting Cell Ranger multi Web Summary Files for Fixed RNA Profiling
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cdn.10xgenomics.com/image/upload/v1706742309/support-documents/CG000729_TechNote_WebSummary_Chromium_FixedRNAProfiling_RevA.pdf


#Beginning the argparse module to extract command-line arguments
parser = argparse.ArgumentParser(description="This is a script that will summarize cellranger count results from at least one cellranger output and create a summary within the 4_reports directory of the project. It accepts a data directory that contains at least one cellranger count results.")
#Creating the following optional arguments; some have default values
parser.add_argument('--dir', type=dir_path, help='Data directory path that contains individually named cellranger count results for samples', required=True)
parser.add_argument('--outdir', type=dir_path, help='Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it.', required=True)
parser.add_argument('--genome', type=str, help='Only specify the genome you want to recover data from a multiple genome alignment', required=False)

#Converts argument strings to objects and assigns them as attributes of the namespace; e.g. --id -> args.id
args = parser.parse_args()


MasterDF = pandas.DataFrame()

for filename in glob.glob(os.path.join(args.dir, "*","metrics_summary_updated.csv")):
    # print(filename)
    
    df = pandas.read_csv(filename)
    
    # Clean the data by removing commas and percent signs
    df = df.replace(",", "", regex=True)  # Remove commas for numbers like '33,635'
    df = df.replace("%", "", regex=True)  # Remove percent signs from values like '81.37%'

    # Convert columns to float (assuming all columns are numeric after cleaning)
    # If only specific columns need conversion, specify the column names like: df["Metric"] = df["Metric"].astype('float')
    df = df.astype('float', errors='ignore')  # Ignore conversion errors for non-numeric columns

    
    SampleID = filename.split("/outs/per_sample_outs/")[1].split("/")[0]
    print(SampleID)
    df["Sample ID"] = SampleID
    
    
    Warnings = ""
    MajorWarnings = ""
    TotalWarnings = 0
    
    
    # Assuming you already have a DataFrame 'df'
    # Clean the "Cells" column by removing any non-numeric characters (e.g., commas)
    df["Cells"] = df["Cells"].replace({",": ""}, regex=True)  # Remove commas
    
    # Convert the "Cells" column to a numeric type (int or float)
    df["Cells"] = pd.to_numeric(df["Cells"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN

    if args.genome == None:
        if df.iloc[0]["Cells"] < 500:
            Warnings = Warnings + "Cells < 500, "
            TotalWarnings += 1
    elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Estimated Number of Cell Partitions"] < 500:
            Warnings = Warnings + "GRCh38 Estimated Number of Cell Partitions < 500, "
            TotalWarnings += 1
    elif args.genome == "mm10":
        if df.iloc[0]["mm10 Estimated Number of Cell Partitions"] < 500:
            Warnings = Warnings + "mm10 Estimated Number of Cell Partitions < 500, "
            TotalWarnings += 1

    
    if args.genome == None:
        if df.iloc[0]["Cells"] < 100:
            MajorWarnings = MajorWarnings + "Cells < 100, "
            TotalWarnings += 1
    elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Estimated Number of Cell Partitions"] < 100:
            MajorWarnings = MajorWarnings + "GRCh38 Estimated Number of Cell Partitions < 100, "
            TotalWarnings += 1
    elif args.genome == "mm10":
        if df.iloc[0]["mm10 Cells"] < 100:
            MajorWarnings = MajorWarnings + "mm10 Cells < 100, "
            TotalWarnings += 1
    
    
    df["Mean reads per cell"] = df["Mean reads per cell"].replace({",": ""}, regex=True)  # Remove commas
    df["Mean reads per cell"] = pd.to_numeric(df["Mean reads per cell"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN


    if df.iloc[0]["Mean reads per cell"] < 10000:
        MajorWarnings = MajorWarnings + "Mean reads per cell < 10000, "
        TotalWarnings += 1
    
       
    df["Valid barcodes"] = df["Valid barcodes"].replace({",": ""}, regex=True)  # Remove commas
    df["Valid barcodes"] = pd.to_numeric(df["Valid barcodes"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
   

    if df.iloc[0]["Valid barcodes"] < 75:
        Warnings = Warnings + "Valid barcodes < 75%, "
        TotalWarnings += 1
    
    df["Q30 barcodes"] = df["Q30 barcodes"].replace({",": ""}, regex=True)  # Remove commas
    df["Q30 barcodes"] = pd.to_numeric(df["Q30 barcodes"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
   
        
    if df.iloc[0]["Q30 barcodes"] < 55:
        Warnings = Warnings + "Q30 barcodes < 55%, "
        TotalWarnings += 1
  
    df["Q30 RNA read"] = df["Q30 RNA read"].replace({",": ""}, regex=True)  # Remove commas
    df["Q30 RNA read"] = pd.to_numeric(df["Q30 RNA read"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
     
    if df.iloc[0]["Q30 RNA read"] < 65:
        Warnings = Warnings + "Q30 RNA read < 65%, "
        TotalWarnings += 1
    
      
    df["Q30 UMI"] = df["Q30 UMI"].replace({",": ""}, regex=True)  # Remove commas
    df["Q30 UMI"] = pd.to_numeric(df["Q30 UMI"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
      
    if df.iloc[0]["Q30 UMI"] < 75:
        Warnings = Warnings + "Q30 UMI < 75%, "
        TotalWarnings += 1
    
          
    df["Reads confidently mapped to probe set"] = df["Reads confidently mapped to probe set"].replace({",": ""}, regex=True)  # Remove commas
    df["Reads confidently mapped to probe set"] = pd.to_numeric(df["Reads confidently mapped to probe set"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
      
    if args.genome == None:
        if df.iloc[0]["Reads confidently mapped to probe set"] < 50:
            Warnings = Warnings + "Reads confidently mapped to probe set < 50%, "
            TotalWarnings += 1
    elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Reads confidently mapped to probe set"] < 50:
            Warnings = Warnings + "GRCh38 Reads confidently mapped to probe set < 50%, "
            TotalWarnings += 1
    elif args.genome == "mm10":
        if df.iloc[0]["mm10 Reads confidently mapped to probe set"] < 50:
            Warnings = Warnings + "mm10 Reads confidently mapped to probe set < 50%, "
            TotalWarnings += 1
            
            
         
    df["Reads confidently mapped to filtered probe set"] = df["Reads confidently mapped to filtered probe set"].replace({",": ""}, regex=True)  # Remove commas
    df["Reads confidently mapped to filtered probe set"] = pd.to_numeric(df["Reads confidently mapped to filtered probe set"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
      
    if args.genome == None:
        if df.iloc[0]["Reads confidently mapped to filtered probe set"] < 50:
            MajorWarnings = MajorWarnings + "Reads confidently mapped to filtered probe set < 50%, "
            TotalWarnings += 1
    elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Reads confidently mapped to filtered probe set"] < 50:
            MajorWarnings = MajorWarnings + "GRCh38 Reads confidently mapped to filtered probe set < 50%, "
            TotalWarnings += 1

    elif args.genome == "mm10":
        if df.iloc[0]["mm10 Reads confidently mapped to filtered probe set"] < 50:
            MajorWarnings = MajorWarnings + "mm10 Reads confidently mapped to filtered probe set < 50%, "
            TotalWarnings += 1

         
    df["Reads mapped to probe set"] = df["Reads mapped to probe set"].replace({",": ""}, regex=True)  # Remove commas
    df["Reads mapped to probe set"] = pd.to_numeric(df["Reads mapped to probe set"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN

    if args.genome == None:
        if df.iloc[0]["Reads mapped to probe set"] < 50:
            MajorWarnings = MajorWarnings + "Reads mapped to probe set < 50%, "
            TotalWarnings += 1
    elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Reads mapped to probe set"] < 50:
            MajorWarnings = MajorWarnings + "GRCh38 Reads mapped to probe set < 50%, "
            TotalWarnings += 1

    elif args.genome == "mm10":
        if df.iloc[0]["mm10 Reads mapped to probe set"] < 50:
            MajorWarnings = MajorWarnings + "mm10 Reads mapped to probe set < 50%, "
            TotalWarnings += 1



    df["Fraction of initial cell barcodes passing high occupancy GEM filtering"] = df["Fraction of initial cell barcodes passing high occupancy GEM filtering"].replace({",": ""}, regex=True)  # Remove commas
    df["Fraction of initial cell barcodes passing high occupancy GEM filtering"] = pd.to_numeric(df["Fraction of initial cell barcodes passing high occupancy GEM filtering"], errors="coerce")  # This will convert strings to numbers, and non-convertible strings will be NaN
   

    if args.genome == None:
        if df.iloc[0]["Fraction of initial cell barcodes passing high occupancy GEM filtering"] < 90:
            MajorWarnings = MajorWarnings + "Fraction of initial cell barcodes passing high occupancy GEM filtering < 90%, "
            TotalWarnings += 1
    elif args.genome == "GRCh38":
        if df.iloc[0]["GRCh38 Fraction of initial cell barcodes passing high occupancy GEM filtering"] < 90:
            MajorWarnings = MajorWarnings + "GRCh38 Fraction of initial cell barcodes passing high occupancy GEM filtering < 90%, "
            TotalWarnings += 1
    elif args.genome == "mm10":
        if df.iloc[0]["mm10 Fraction of initial cell barcodes passing high occupancy GEM filtering"] < 90:
            MajorWarnings = MajorWarnings + "mm10 Fraction of initial cell barcodes passing high occupancy GEM filtering < 90%, "
            TotalWarnings += 1
 
            
    df["Warnings"] = Warnings
    df["MajorWarnings"] = MajorWarnings
    df["Total Warnings"] = TotalWarnings
        
    MasterDF = pandas.concat([MasterDF, df])

MasterDF.to_csv( args.outdir + "QC_Summary_CellRanger_Report.tsv", sep = "\t", index = False)
