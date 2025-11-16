# Databricks notebook source
# MAGIC %md
# MAGIC ## [Create](https://docs.databricks.com/en/sql/language-manual/sql-ref-syntax-ddl-create-sql-function.html) and register `scientific2simple` as a [`UDF`](https://docs.databricks.com/en/udf/unity-catalog.html) to Unity Catalog 
# MAGIC To democratize the usage or our [Classified Proteins]($./02_TransformerCNN_Protein_Classification) in our [AI/BI Genie](https://www.databricks.com/product/ai-bi/genie) [Space](https://docs.databricks.com/en/genie/index.html), we can leverage [Foundational Models API](https://docs.databricks.com/en/machine-learning/foundation-models/index.html#use-foundation-apis) to help simplify _protein scientific terms_ and provide _layman simple terms with definition_.
# MAGIC
# MAGIC To do so, we will use [`ai_query()`](https://docs.databricks.com/en/sql/language-manual/functions/ai_query.html) [to query the Foundation Model](https://docs.databricks.com/en/large-language-models/how-to-ai-query.html) `databricks-meta-llama-3-3-70b-instruct`  in a specific way to bulk convert the `OrganismName` scientific terms to layman simple terms and to provide [`UDF`](https://docs.databricks.com/en/udf/unity-catalog.html) [structured outputs](https://docs.databricks.com/en/machine-learning/model-serving/structured-outputs.html). 
# MAGIC
# MAGIC To reuse the logic of converting scientific terms to simpler layman terms (e.g. as a tool in [function-calling](https://docs.databricks.com/en/machine-learning/model-serving/function-calling.html)), we can register the the `ai_query()` as a [`UDF function`](https://docs.databricks.com/en/udf/unity-catalog.html) and then use it in a query to create the desired table with simplified terms and their meanings.
# MAGIC
# MAGIC We'll walk through how this can be achieved.
# MAGIC
# MAGIC _NB: we can use a serverless compute for setting this up._ 
# MAGIC
# MAGIC <!-- 
# MAGIC https://docs.databricks.com/en/genie/trusted-assets.html
# MAGIC
# MAGIC https://docs.databricks.com/en/sql/language-manual/sql-ref-syntax-ddl-create-function.html
# MAGIC
# MAGIC https://docs.databricks.com/en/udf/unity-catalog.html
# MAGIC  -->

# COMMAND ----------

# DBTITLE 1,run notebook utils
# MAGIC %run ./utils

# COMMAND ----------

# MAGIC %md
# MAGIC ### [1] Let's review our classified proteins data: `proteinclassification_tiny` 
# MAGIC

# COMMAND ----------

# DBTITLE 1,UC variables
remove_widgets() 
uc_config = setup_uc_paths(spark=None, use_widgets=False); ## if you update the values in widgets -- it will automatically trigger an update of the UC paths

# Extract catalog, schema, volume names
catalog_name = uc_config["catalog_name"]
schema_name = uc_config["schema_name"]
volume_name = uc_config["volume_name"]

# COMMAND ----------

# DBTITLE 1,Read in proteinclassification_tiny
sDF = spark.table(f"{catalog_name}.{schema_name}.proteinclassification_tiny")

# COMMAND ----------

# DBTITLE 1,View sparkDataFrame
display(sDF.sort("OrganismName"))

# display(sDF.groupby('OrganismName', 'ProteinName', 'ID').count().sort('OrganismName'))

# COMMAND ----------

# MAGIC %md
# MAGIC We can see that some Organisms are associated with different types of (e.g soluble / membrane transport) proteins as denoted by their corresponding `ProteinNames`, `OrganismIdentifier`, `GeneName`, and `IDs`

# COMMAND ----------

# DBTITLE 1,OrganismName -- Scientific Geneology
from pyspark.sql import functions as F, types as T

display(sDF.groupby('OrganismName').agg((F.size(F.collect_set('ProteinName')).alias('N_uniqueProteins')),
                                        F.collect_set('ProteinName'),                                          
                                        F.collect_set('ID')
                                       ).sort('N_uniqueProteins', ascending=False))

# COMMAND ----------

# MAGIC %md
# MAGIC It would be helpful if we could convert these scientific `OrganismName` into simple layman terms so that non-SME biologists can more easily explore the data. 
# MAGIC
# MAGIC We will write out the unique set of `OrganismName` as a separate Delta Table to our Unity Catalog to use for exploring the conversion of scientific to simple terms using LLMs. The derived simple terms can then be linked by by the scientific `OrganismName` later. 

# COMMAND ----------

# DBTITLE 1,Write/Read tinysample_organism_info
# sDF.select('OrganismName').distinct().write.mode("overwrite").option("mergeSchema", "true").saveAsTable(f"{catalog_name}.{schema_name}.tinysample_organism_info")

## check/read tinysample_organism_info 
sDF_orginfo = spark.table(f"{catalog_name}.{schema_name}.tinysample_organism_info")
display(sDF_orginfo)

# COMMAND ----------

# MAGIC %md
# MAGIC ### [2] Use `ai_query()` to leverage LLM to help simplify scientific terms
# MAGIC
# MAGIC We can query and look up these scientific `OrganismName` with the help of LLMs and try to organize the query output in a way that will facilitate our use of the corresponding protein dataset. 
# MAGIC
# MAGIC Let's use the [Foundation Model](https://docs.databricks.com/en/large-language-models/how-to-ai-query.html) `databricks-meta-llama-3-1-70b-instruct` model and guide it to help our case, by providing a detailed prompt. We will wrap it in an [`ai_query()`](https://docs.databricks.com/en/sql/language-manual/functions/ai_query.html) that request that it provides a [structured output](https://docs.databricks.com/en/machine-learning/model-serving/structured-outputs.html) JSON response:
# MAGIC
# MAGIC
# MAGIC e.g. 
# MAGIC
# MAGIC | Scientific `OrganismName` | `SimpleTermDict` (JSON Response) | 
# MAGIC |-----------------|---------------|
# MAGIC | Danio rerio     | {"simple_term": "Zebrafish", "meaning": "A species of freshwater fish used as a model organism in scientific research, particularly in the fields of developmental biology and neurology."} | 
# MAGIC
# MAGIC
# MAGIC <!-- | Scientific `OrganismName` | `SimpleTermDict` (JSON Response) | `simple_term` | `meaning` |
# MAGIC |-----------------|---------------|-------------|---------|
# MAGIC | Danio rerio     | {"simple_term": "Zebrafish", "meaning": "A species of freshwater fish used as a model organism in scientific research, particularly in the fields of developmental biology and neurology."} | Zebrafish | A species of freshwater fish used as a model organism in scientific research, particularly in the fields of developmental biology and neurology. | -->

# COMMAND ----------

# DBTITLE 1,TEST: Can LLM help simplify scientific terms?
# MAGIC %sql
# MAGIC -- Use the ai_query function directly in a query and perform the necessary scientific2simple term transformations in a single query
# MAGIC
# MAGIC -- CREATE OR REPLACE TABLE ${catalog_name}.${schema_name}.tinysample_organism_info_scientificNsimple USING DELTA AS
# MAGIC SELECT
# MAGIC   OrganismName,
# MAGIC   SimpleTermDict,
# MAGIC   get_json_object(SimpleTermDict, '$.simple_term') AS Organism_SimpleTerm,
# MAGIC   get_json_object(SimpleTermDict, '$.meaning') AS Organism_Definition
# MAGIC FROM (
# MAGIC   SELECT
# MAGIC     OrganismName,
# MAGIC     ai_query(
# MAGIC       'databricks-meta-llama-3-3-70b-instruct', 
# MAGIC       -- please update to model and/or version that exists in your workspace e.g.
# MAGIC       -- 'databricks-claude-sonnet-4-5',
# MAGIC       CONCAT(
# MAGIC         'As a knowledgeable and factual encyclopedia you respond succinctly. Provide a dictionary of responses in the format {"key": "response"} to the following keys a "simple layman term" and "meaning" for the given "scientific term": ', 
# MAGIC         OrganismName, 
# MAGIC         '. Output just the dictionary {"simple_term": "response", "meaning": "response"}.'
# MAGIC       )
# MAGIC     ) AS SimpleTermDict
# MAGIC   FROM mmt_demos2.ai_driven_drug_discovery.tinysample_organism_info --full name required
# MAGIC ) ORDER BY OrganismName

# COMMAND ----------

# MAGIC %md
# MAGIC ### [3] `ai_query()` structured output JSON parsed: 
# MAGIC
# MAGIC Our [`ai_query()`](https://docs.databricks.com/en/sql/language-manual/functions/ai_query.html) outputs the keys `simple_term` and `meaning` of the scientific `OrganismName` and their corresponding values, which we extract using `get_json_object` parsing within the single SQL query. 
# MAGIC
# MAGIC **To reuse the logic of converting scientific terms to simpler layman terms (e.g. as a tool in [function-calling](https://docs.databricks.com/en/machine-learning/model-serving/function-calling.html)), we can**    
# MAGIC **-A. register the the `ai_query()` as a [`UDF function`](https://docs.databricks.com/en/udf/unity-catalog.html) and then**    
# MAGIC **-B. use it in a query to create the desired table with simplified terms and their meanings.**  
# MAGIC
# MAGIC (_This is particularly relevant given that we cannot directly include `get_json_object` parsing within the `SQL function` definition itself, we can create a view or a separate query that uses the UC registered `sql function` and then apply `get_json_object` to parse the JSON response._)
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC #### [3A] Define an `ai_query()` to process `scientific2simple_dict` as `UDF function` + register to UC
# MAGIC
# MAGIC You can use SQL and pyspark.sql to define and register the UDF function to Unity Catalog.  

# COMMAND ----------

# DBTITLE 1,[A] Define & Register UDF function
### via PySpark 
# Define the SQL function
create_function_query = f"""
CREATE OR REPLACE FUNCTION {catalog_name}.{schema_name}.scientific2simple_dict(OrganismName STRING)
RETURNS STRING
COMMENT 'Returns a dictionary of simple terms and definitions for the given scientific term as JSON output.'
LANGUAGE SQL
RETURN ai_query(
  'databricks-meta-llama-3-3-70b-instruct',
  CONCAT('As a knowledgeable and factual encyclopedia you respond succinctly. Provide a dictionary of responses in the format {{"key": "response"}} to the following keys a "simple layman term" and "meaning" for the given "scientific term": ', OrganismName, '. Output just the dictionary {{"simple_term": "response", "meaning": "response"}}.')
);
"""

# Execute the query to create the function
spark.sql(create_function_query)


#-------------------------------------------------------------------------------------------------------------------------------

### via SQL 
# %sql

# CREATE OR REPLACE FUNCTION ${catalog_name}.${schema_name}.scientific2simple_dict(OrganismName STRING)
# RETURNS STRING
# COMMENT 'Returns a dictionary of simple terms and definitions for the given scientific term as JSON output.'
# LANGUAGE SQL
# RETURN ai_query(
#   'databricks-meta-llama-3-3-70b-instruct',
#   CONCAT('As a knowledgeable and factual encyclopedia you respond succinctly. Provide a dictionary of responses in the format {"key": "response"} to the following keys a "simple layman term" and "meaning" for the given "scientific term": ', OrganismName, '. Output just the dictionary {"simple_term": "response", "meaning": "response"}.')
# );

### -- Test the registered function
# SELECT
#   OrganismName,
#   ${catalog_name}.${schema_name}.scientific2simple_dict(OrganismName) AS SimpleTermDict
# FROM ${catalog_name}.${schema_name}.tinysample_organism_info;


# COMMAND ----------

# MAGIC %md
# MAGIC #### [3B] Use the registered SQL Function in a Query 
# MAGIC Additionally, perform necessary transformations to extract the simplified scientific `OrganismName` and definitions, then write out transformed data to UC. 
# MAGIC
# MAGIC NB: You can definitely leverage available LLMs to simply extract everything (human layman terms and any other downstream probe on the data) at once in an interactive query in genie but if these are common queries on a bulk of data it would be most efficient and cost-effective to preprocess the joins prior. 

# COMMAND ----------

# DBTITLE 1,[B] Use SQL function via Spark SQL for Batch Inferencing
# Use a Common Table Expression (CTE) or a subquery to ensure that the ai_query function is called only once per OrganismName

orginfo_sDF = spark.sql(f"""
WITH SimpleTermData AS (
  SELECT
    OrganismName,
    {catalog_name}.{schema_name}.scientific2simple_dict(OrganismName) AS SimpleTermDict
  FROM {catalog_name}.{schema_name}.tinysample_organism_info
)
SELECT
  OrganismName,
  SimpleTermDict,
  get_json_object(SimpleTermDict, '$.simple_term') AS Organism_SimpleTerm,
  get_json_object(SimpleTermDict, '$.meaning') AS Organism_Definition
FROM SimpleTermData
-- ORDER BY OrganismName; -- this slows down the query
""")

# writing out the query to a UC Delta table speeds up subsequent queries and can be used for incremental updates
orginfo_sDF.write.mode("overwrite").option("mergeSchema", "true").saveAsTable(f"{catalog_name}.{schema_name}.tinysample_organism_info_scientificNsimple")

# COMMAND ----------

# MAGIC %md
# MAGIC ### [4] Read and display the resulting saved table with simplified `OrganismName`

# COMMAND ----------

# DBTITLE 1,Resulting simplified OrganismName
orginfo_sDF = spark.table(f"{catalog_name}.{schema_name}.tinysample_organism_info_scientificNsimple")

display(orginfo_sDF)
