# Databricks notebook source
# DBTITLE 1,Install necessary packages.
# MAGIC %pip install -q biopython
# MAGIC dbutils.library.restartPython()
# MAGIC
# MAGIC ## For Spark Declarative Pipelines, dependecies need to be installed right at the top of the notebook

# COMMAND ----------

# MAGIC %md
# MAGIC # Protein Data Pre-Processing using Spark Declarative Pipelines

# COMMAND ----------

# MAGIC %md
# MAGIC For our [AI/BI Dashboard](https://docs.databricks.com/en/dashboards/index.html) with [enabled Genie Space](https://docs.databricks.com/en/dashboards/index.html#enable-a-genie-space-from-your-dashboard) example we will use the [UniProt FASTA protein sequence data](https://www.uniprot.org/help/downloads), which we downloaded to the [Unity Catalog](https://docs.databricks.com/en/catalogs/index.html) [Volumes](https://docs.databricks.com/en/volumes/index.html) associated with our project [Schema](https://docs.databricks.com/en/schemas/index.html#what-is-a-schema): `<catalog_name>.ai_driven_drug_discovery`
# MAGIC As with most `raw` data, we will need to pre-process to 'clean' and extract relevant information before we can interact and use the data for our downstream explorations. 
# MAGIC
# MAGIC Here we will define a `Protein Data Processing` [Spark Declarative Pipelines](https://www.databricks.com/product/delta-live-tables) to help extract the protein information e.g. *`identifiers`*, *`sequences`*, *`names`*, *`organism information`*, *`gene names`*, *`protein existence`*, and perform any additional transformation(s) e.g. calculating *`molecular weight`* in a sequence of tasks.   
# MAGIC         

# COMMAND ----------

# MAGIC %md        
# MAGIC
# MAGIC  <img src="../assets/imgs/ProteinDataProcessing_DLT.png" alt="Protein Data Processing DLT" />

# COMMAND ----------

# DBTITLE 1,in case of re-establishing the DLT
# # List of tables to drop
# tables_to_drop = [
#     f"{catalog_name}.{schema_name}.bronze_protein",
#     f"{catalog_name}.{schema_name}.silver_protein",
#     f"{catalog_name}.{schema_name}..enriched_protein"
# ]

# # Drop each table
# for table in tables_to_drop:
#     spark.sql(f"DROP TABLE IF EXISTS {table}")

# COMMAND ----------

# MAGIC %md
# MAGIC

# COMMAND ----------

# DBTITLE 1,Define the DLT pipeline
# MAGIC %md
# MAGIC ## Define Tasks for our DLT pipeline
# MAGIC
# MAGIC <!-- [Protein Data Processing](https://e2-demo-west.cloud.databricks.com/pipelines/c6a3e57b-1c44-476f-9e56-c8e05d9975f5/updates/2adf25f7-deb6-40fa-8f87-68efd0ca05a0?o=2556758628403379%3Fparent%3Dfolders%2F1625373258638091)  -->

# COMMAND ----------

# DBTITLE 1,Import required packages/libraries.
import dlt
from Bio import SeqIO
from pyspark.sql import Row
import io

# Additional library/package import will be included in cells related to the DLT pipeline Tasks specified below

# COMMAND ----------

# MAGIC %md
# MAGIC ### [1] Create `bronze_protein` DLT 
# MAGIC Load 500,000 protein sequences from 1 text file into Bronze Table

# COMMAND ----------

# DBTITLE 1,Define UC variables
# Use widgets for dynamic configuration
dbutils.widgets.text("catalog_name", "mmt_demos2", "Catalog Name")
dbutils.widgets.text("schema_name", "ai_driven_drug_discovery", "Schema Name")
dbutils.widgets.text("volume_name", "protein_seq", "Volume Name")

catalog_name = dbutils.widgets.get("catalog_name")
schema_name = dbutils.widgets.get("schema_name")
volume_name = dbutils.widgets.get("volume_name")

print(f"Using catalog: {catalog_name}, schema: {schema_name}, volume: {volume_name}")

# COMMAND ----------

# DBTITLE 1,Read in FASTA file as DLT materialize table
# Initialize lists to hold the data
records = []

# Read from UC Volumes and Parse the FASTA file

file_path = f'/Volumes/{catalog_name}/{schema_name}/{volume_name}/uniprot_sprot.fasta'

with open(file_path, "r") as f:
     for record in SeqIO.parse(f, "fasta"):
      id = record.id
      sequence = str(record.seq)
      description = record.description
      
      records.append(Row(ID=id, Sequence=sequence, Description=description))

df = spark.createDataFrame(records)


##Create a Bronze Delta Live Table
@dlt.expect_or_drop("Empty ID","ID != '' ")
@dlt.table(
    comment = "FASTA Data",
    table_properties = {"quality": "bronze"})    
def bronze_protein():
    return spark.createDataFrame(records)


# COMMAND ----------

# MAGIC %md
# MAGIC ### [2] Create `silver_protein` DLT 
# MAGIC Parse through the text and format data into columns

# COMMAND ----------

# DBTITLE 1,Extract Protein Info. from FASTA file
## Import require functions 
from pyspark.sql.functions import regexp_extract

## Create a Silver Delta Live Table 
@dlt.expect_or_drop("Empty ID","ID != '' ")
@dlt.table
def silver_protein():
    # Regular expressions for each field
    os_regex = r'OS=([^ ]+ [^ ]+|\([^)]+\))'
    ox_regex = r'OX=(\d+)'
    gn_regex = r'GN=([^ ]+)'
    pe_regex = r'PE=(\d)'
    sv_regex = r'SV=(\d)'


    fasta_df = dlt.read("bronze_protein")

    # Extract ProteinName
    fasta_df = fasta_df.withColumn("ProteinName", regexp_extract("Description", r" (.+?) OS=", 1))

    # Extract and create new columns for OrganismName, OrganismIdentifier, GeneName, ProteinExistence, SequenceVersion
    fasta_df = fasta_df.withColumn('OrganismName', regexp_extract('Description', os_regex, 1))
    fasta_df = fasta_df.withColumn('OrganismIdentifier', regexp_extract('Description', ox_regex, 1))
    fasta_df = fasta_df.withColumn('GeneName', regexp_extract('Description', gn_regex, 1))
    fasta_df = fasta_df.withColumn('ProteinExistence', regexp_extract('Description', pe_regex, 1))
    fasta_df = fasta_df.withColumn('SequenceVersion', regexp_extract('Description', sv_regex, 1))

    return fasta_df


# COMMAND ----------

# MAGIC %md
# MAGIC ### [3] Create `enriched_protein` DLT
# MAGIC  We can also demonstrate how fast and easy it is to use third party libraries e.g. ```Bio.SeqUtils' molecular_weight``` within a [`Pandas` User-Defined Function](https://docs.databricks.com/en/udf/pandas.html) inside a DLT pipeline task to calculate ***molecular weights*** of each molecule in a vectorized process.    
# MAGIC       
# MAGIC  Here, we are calculating *molecular weights* for 500,000 molecules and this takes about 30 seconds using a serverless DLT

# COMMAND ----------

# DBTITLE 1,Include Molecular Weight Calculation
from pyspark.sql.functions import pandas_udf
from pyspark.sql.types import DoubleType
from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq
import pandas as pd

# Assume Arrow is already enabled in your Databricks cluster

# Define our Pandas UDF to calculate molecular weight using Bio.SeqUtils' molecular_weight module
@pandas_udf(DoubleType())
def get_molecular_weight_pandas_udf(sequence: pd.Series) -> pd.Series:
    def calculate_mw(seq):
        try:
            # Attempt to calculate the molecular weight
            return molecular_weight(Seq(seq), seq_type="protein")
        except ValueError as e:
            return 1.0
                
    
    return sequence.apply(calculate_mw)

# DLT Pipeline function to enrich and modify the table
@dlt.expect_or_drop("Empty ID","ID != '' ")
@dlt.table
def enriched_protein():
    # Load the existing silver_protein_parsed table
    df = dlt.read("silver_protein")
    
    # Add the "Molecular Weight" column using the Pandas UDF to vectorize the molecular weights calculation
    df = df.withColumn("Molecular_Weight", get_molecular_weight_pandas_udf(df["sequence"]))
    
    # Drop the "Description" column from the DataFrame
    df = df.drop("Description")
    
    return df


# COMMAND ----------

# MAGIC %md
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC ### [4] Review our defined DLT pipeline & `enriched_protein` sparkDF
# MAGIC - Our defined DLT pipeline: [Protein Data Processing](https://e2-demo-field-eng.cloud.databricks.com/pipelines/f7cfd697-0adf-4719-85d0-26779034145b/updates/14a457fe-215b-4f6b-99d1-c5d81ce3345d?o=1444828305810485)    
# MAGIC
# MAGIC - Use ``` (Ctrl + `)``` to view Terminal, DLT, Logs within notebook  
# MAGIC - Turn on ETL Pipeline editor 
