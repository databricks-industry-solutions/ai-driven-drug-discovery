# Databricks notebook source
# DBTITLE 1,Install dependencies
!pip install biopython
# !pip install s3fs

# COMMAND ----------

# DBTITLE 1,Load libraries
import dlt
from Bio import SeqIO
from pyspark.sql import Row
import io
# import s3fs

# Additional library/package import will be included in cells related to the DLT pipeline Tasks specified below

# COMMAND ----------

# DBTITLE 1,run utils
# MAGIC %run ./utils

# COMMAND ----------

# DBTITLE 1,define UC variables in utils and extract them
remove_widgets() 
uc_config = setup_uc_paths(spark=None, use_widgets=True); ## if you update the values in widgest -- it will automatically trigger an update of the UC paths

# Extract catalog, schema, volume names
catalog_name = uc_config["catalog_name"]
schema_name = uc_config["schema_name"]
volume_name = uc_config["volume_name

# COMMAND ----------

# DBTITLE 1,Define UC variables
# catalog_name = "mmt_demos2"
# schema_name = "ai_driven_drug_discovery"


# # CATALOG_NAME = "<catalog_name>"              # TODO: Replace with your Unity Catalog name
# # SCHEMA_NAME = "ai_driven_drug_discovery"     # TODO: Replace with your schema name
# # VOLUME_NAME = "protein_seq"                  # TODO: Replace with your volume name for file storage
# # ENDPOINT_NAME = "az_openai_gpt4o"            # TODO: Replace with your AI Gateway endpoint name


# COMMAND ----------

# DBTITLE 1,Read in FASTA file as DLT materialize table
from pyspark.sql import Row
from Bio import SeqIO

# Initialize lists to hold the data
records = []

# Read from DBFS and Parse the FASTA file
file_path = f'/Volumes/{catalog_name}/{schema_name}/protein_seq/uniprot_sprot.fasta'

with open(file_path, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        id = record.id
        sequence = str(record.seq)
        description = record.description
        records.append(Row(ID=id, Sequence=sequence, Description=description))

df = spark.createDataFrame(records)

# Save as Bronze Delta Table
df.write.format("delta").mode("overwrite").saveAsTable(f"{catalog_name}.{schema_name}.bronze_protein")

# COMMAND ----------

# DBTITLE 1,Extract Protein Info. from FASTA file
from pyspark.sql.functions import regexp_extract

# Load the existing bronze_protein table
fasta_df = spark.table(f"{catalog_name}.{schema_name}.bronze_protein")

# Regular expressions for each field
os_regex = r'OS=([^ ]+ [^ ]+|\([^)]+\))'
ox_regex = r'OX=(\d+)'
gn_regex = r'GN=([^ ]+)'
pe_regex = r'PE=(\d)'
sv_regex = r'SV=(\d)'

# Extract ProteinName
fasta_df = fasta_df.withColumn("ProteinName", regexp_extract("Description", r" (.+?) OS=", 1))

# Extract and create new columns for OrganismName, OrganismIdentifier, GeneName, ProteinExistence, SequenceVersion
fasta_df = fasta_df.withColumn('OrganismName', regexp_extract('Description', os_regex, 1))
fasta_df = fasta_df.withColumn('OrganismIdentifier', regexp_extract('Description', ox_regex, 1))
fasta_df = fasta_df.withColumn('GeneName', regexp_extract('Description', gn_regex, 1))
fasta_df = fasta_df.withColumn('ProteinExistence', regexp_extract('Description', pe_regex, 1))
fasta_df = fasta_df.withColumn('SequenceVersion', regexp_extract('Description', sv_regex, 1))

# Save as Silver Delta Table
fasta_df.write.format("delta").mode("overwrite").saveAsTable(f"{catalog_name}.{schema_name}.silver_protein")

# COMMAND ----------

# DBTITLE 1,Include Molecular Weight Calculation
from pyspark.sql.functions import pandas_udf
from pyspark.sql.types import DoubleType
from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq
import pandas as pd

# Define our Pandas UDF to calculate molecular weight using Bio.SeqUtils' molecular_weight module
@pandas_udf(DoubleType())
def get_molecular_weight_pandas_udf(sequence: pd.Series) -> pd.Series:
    def calculate_mw(seq):
        try:
            return molecular_weight(Seq(seq), seq_type="protein")
        except ValueError as e:
            return 1.0
    return sequence.apply(calculate_mw)

# Load the existing silver_protein table
df = spark.table(f"{catalog_name}.{schema_name}.silver_protein")

# Add the "Molecular Weight" column using the Pandas UDF to vectorize the molecular weights calculation
df = df.withColumn("Molecular_Weight", get_molecular_weight_pandas_udf(df["Sequence"]))

# Drop the "Description" column from the DataFrame
df = df.drop("Description")

# Save as Enriched Delta Table
df.write.format("delta").mode("overwrite").saveAsTable(f"{catalog_name}.{schema_name}.enriched_protein")

# COMMAND ----------

read_enriched_protein = spark.table(f"{catalog_name}.{schema_name}.enriched_protein")

display(read_enriched_protein)
