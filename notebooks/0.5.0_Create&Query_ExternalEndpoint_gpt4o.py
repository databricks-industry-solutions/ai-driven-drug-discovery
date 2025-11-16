# Databricks notebook source
# MAGIC %md
# MAGIC ## Create a [Mosaic AI Model Serving](https://learn.microsoft.com/en-us/azure/databricks/machine-learning/model-serving/) endpoint to serve an [`External Model`](https://docs.databricks.com/en/generative-ai/external-models/index.html)  
# MAGIC
# MAGIC [`At the time of this Solution's development, GPT-4o just surfaced. Today you have many other Foundational Models to consider. However, we use Azure GPT-4o here to illustrated how to create an external model serving endpoint.`]
# MAGIC
# MAGIC In order to leverage `GPT-4o's Medical and Scientific knowledge capabilities` (ref to `Health` and `Scientific Capabilities` sections under `Social Impacts` in [`gpt-4o-system-card`](https://openai.com/index/gpt-4o-system-card/)), we will serve the external `Azure Openai GPT-4o` model endpoint to use with our [AIBI Genie Space](https://docs.databricks.com/en/genie/index.html) 
# MAGIC
# MAGIC (_Omit steps `1-3` if your workspace already has an external endpoint serving `gpt-4o` with which you can call and make inference_)

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC There are different ways to create external endpoints to serve `{Azure}OpenAI/non-Databricks` models e.g. via the [UI](https://docs.databricks.com/en/machine-learning/model-serving/create-foundation-model-endpoints.html#ext-model-endpoint) or[`mlflow.deployments`](https://docs.databricks.com/en/generative-ai/tutorials/external-models-tutorial.html) (_NB: at the time of documenting, `mlflow.deployments` does not support Mosaic AI Gateway options_).
# MAGIC
# MAGIC Here, we use the [serving-endpoints API](https://docs.databricks.com/en/machine-learning/model-serving/create-foundation-model-endpoints.html#language-REST%C2%A0API) ([ref](https://docs.databricks.com/api/workspace/servingendpoints/create)) to serve the `Azure OpenAI GPT4o` model and include the [Mosaic AI Gateway](https://docs.databricks.com/en/ai-gateway/index.html) options e.g. 
# MAGIC - track endpoint usage and associated costs using [system tables](https://docs.databricks.com/en/admin/system-tables/index.html) 
# MAGIC - include payload logging for model inference data audit using [inference table](https://docs.databricks.com/en/machine-learning/model-serving/inference-tables.html#what-is) which could be further used for monitoring if desired   
# MAGIC
# MAGIC [These options can be configured using the UI or specified via code](https://docs.databricks.com/en/ai-gateway/configure-ai-gateway-endpoints.html).    
# MAGIC
# MAGIC
# MAGIC <!-- # ref 
# MAGIC # https://docs.databricks.com/en/generative-ai/tutorials/external-models-tutorial.html 
# MAGIC
# MAGIC # https://docs.databricks.com/en/machine-learning/model-serving/create-foundation-model-endpoints.html#language-REST%C2%A0API
# MAGIC # https://docs.databricks.com/api/workspace/servingendpoints/create
# MAGIC -->
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC ### [1] The following pre-requisites are needed: 
# MAGIC
# MAGIC - [Azure OpenAI API subscription](https://portal.azure.com/#create/Microsoft.CognitiveServicesOpenAI) (Keys & Endpoint Info.)
# MAGIC   - API Base = Endpoint URL: `https://<served_ai_foundry_endpoint_name>.openai.azure.com/` 
# MAGIC   <!-- - API Base = Endpoint URL: `https://hls-fedemo-azure-openai.openai.azure.com/`  -->
# MAGIC   - API KEY: _to store as [`Databricks CLI`](https://docs.databricks.com/en/dev-tools/cli/index.html) [`secret` within a `scope`](https://docs.databricks.com/en/security/secrets/index.html)_
# MAGIC
# MAGIC - [Azure AI Foundry/Azure OpenAI Service/Deployments](https://learn.microsoft.com/en-us/azure/ai-services/openai/how-to/create-resource?pivots=web-portal#deploy-a-model) 
# MAGIC   - Deploy required external model
# MAGIC     - Deployment Name (under Deployment Info.): `gpt-4o-2024-11-20` (NB model version = `2024-11-20`)
# MAGIC     - Model API Version (from Endpoint Target URI): `2024-08-01-preview`
# MAGIC
# MAGIC     [_NB: at the time this solution was developed these versions were selected -- newer versions likely replaced these over time -- please choose a version that is currently available for your needs._ NOTE also that model version `YYYY-MM-DD` can differ from Model API Version `YYY-MM-DD-{description}`]
# MAGIC
# MAGIC - Workspace [`Personal Access Token`](https://learn.microsoft.com/en-us/azure/databricks/dev-tools/auth/pat) OR [Azure Entra ID Token-based](https://learn.microsoft.com/en-us/azure/databricks/dev-tools/service-prin-aad-token) [`Service Principal`](https://learn.microsoft.com/en-us/azure/databricks/admin/users-groups/service-principals) (recommended for production workloads) 
# MAGIC
# MAGIC <!-- 
# MAGIC Internal Ref: https://docs.google.com/document/d/1Sj-Bz0vRNi6AGYxDBugkoD4Luscevq86Tqio8uLO9GE/edit?pli=1&tab=t.0 
# MAGIC
# MAGIC ## AZ-OPENAI-subscription
# MAGIC # https://portal.azure.com/#@DataBricksInc.onmicrosoft.com/resource/subscriptions/3f2e4d32-8e8d-46d6-82bc-5bb8d962328b/resourceGroups/fe-shared-amer-001/providers/Microsoft.CognitiveServices/accounts/hls-fedemo-azure-openai/overview
# MAGIC
# MAGIC # https://portal.azure.com/#@DataBricksInc.onmicrosoft.com/resource/subscriptions/3f2e4d32-8e8d-46d6-82bc-5bb8d962328b/resourceGroups/fe-shared-amer-001/providers/Microsoft.CognitiveServices/accounts/hls-fedemo-azure-openai/cskeys
# MAGIC
# MAGIC # KEY1: setup scope + secret 
# MAGIC
# MAGIC # endpoint/base: https://hls-fedemo-azure-openai.openai.azure.com/
# MAGIC
# MAGIC ## DEPLOYMENT
# MAGIC # https://ai.azure.com/resource/deployments/%2Fsubscriptions%2F3f2e4d32-8e8d-46d6-82bc-5bb8d962328b%2FresourceGroups%2Ffe-shared-amer-001%2Fproviders%2FMicrosoft.CognitiveServices%2Faccounts%2Fhls-fedemo-azure-openai%2Fdeployments%2Fgpt-4o-2024-11-20?wsid=/subscriptions/3f2e4d32-8e8d-46d6-82bc-5bb8d962328b/resourceGroups/fe-shared-amer-001/providers/Microsoft.CognitiveServices/accounts/hls-fedemo-azure-openai&tid=9f37a392-f0ae-4280-9796-f1864a10effc
# MAGIC
# MAGIC # target uri: https://hls-fedemo-azure-openai.openai.azure.com/openai/deployments/gpt-4o-2024-11-20/chat/completions?api-version=2024-08-01-preview
# MAGIC
# MAGIC # base: https://hls-fedemo-azure-openai.openai.azure.com/
# MAGIC # deployment name: gpt-4o-2024-11-20   ## model version 2024-11-20
# MAGIC # api version: 2024-08-01-preview
# MAGIC -->
# MAGIC

# COMMAND ----------

# DBTITLE 1,Databricks CLI
# MAGIC %md
# MAGIC ### [2] Create Databricks [Secret](https://docs.databricks.com/en/security/secrets/index.html) Scope + Secret for [Azure OpenAI](https://azure.microsoft.com/en-us/products/ai-services/openai-service) [API](https://portal.azure.com/#create/Microsoft.CognitiveServicesOpenAI) KEY
# MAGIC
# MAGIC The [`Databricks CLI`](https://docs.databricks.com/en/dev-tools/cli/index.html) is only supported for interactive use from the web terminal on x86 compute.    
# MAGIC We recommend using    
# MAGIC i) the [notebook web terminal to invoke CLI commands](https://learn.microsoft.com/en-us/azure/databricks/notebooks/notebook-ui#cli) to [create scopes and store secrets](https://docs.databricks.com/en/security/secrets/index.html#create-a-secret); OR     
# MAGIC 2) the [Databricks Python SDK](https://databricks-sdk-py.readthedocs.io/en/latest/workspace/workspace/secrets.html) if you would like to interface with Databricks APIs.
# MAGIC
# MAGIC
# MAGIC We provide some guiding code for Scope and corresponding Secret setup using Databricks CLI here:   
# MAGIC
# MAGIC Check databricks path:
# MAGIC ```
# MAGIC
# MAGIC >> which databricks
# MAGIC /usr/local/bin/databricks
# MAGIC  
# MAGIC ```   
# MAGIC
# MAGIC Check which databricks version ... which prompts the Databricks CLI installion on notebook terminal: 
# MAGIC ```
# MAGIC
# MAGIC >> databricks --version
# MAGIC Installing the CLI...
# MAGIC Installed Databricks CLI v0.237.0 at /root/bin/databricks.
# MAGIC Databricks CLI v0.237.0
# MAGIC  
# MAGIC ```    
# MAGIC
# MAGIC Once installed we can use the CLI to list existing Scopes: 
# MAGIC ```
# MAGIC
# MAGIC >> databricks secrets list-scopes 
# MAGIC  
# MAGIC ```    
# MAGIC
# MAGIC <!-- We will use `scope_name= hls_fedemo_azure_openai` to create a new `scope` and include a `secret` corresponding to the `scope_key= azopenai_api_key` :   -->
# MAGIC We will use **`scope_name`** `= "<{your_scope_name_prefix}_azure_openai>"` to create a new `scope` and include a `secret` corresponding to the **`scope_key`** `= "<azopenai_api_key>"` : 
# MAGIC
# MAGIC ```
# MAGIC
# MAGIC databricks secrets create-scope {scope_name} 
# MAGIC
# MAGIC databricks secrets put-secret {scope_name} {scope_key}
# MAGIC  
# MAGIC ```
# MAGIC ---     
# MAGIC
# MAGIC For our example: 
# MAGIC ```
# MAGIC
# MAGIC databricks secrets create-scope {your_scope_name_prefix}-azure-openai
# MAGIC
# MAGIC databricks secrets put-secret {your_scope_name_prefix}-azure-openai openai_api_key   
# MAGIC ```  
# MAGIC _(This prompts for `KEY` which we input with the `Azure OpenAI API KEY`)_
# MAGIC
# MAGIC
# MAGIC We can check the defined scope and secret using `databricks secrets get-secret {scope}, {key}`.   
# MAGIC (_NB: this outputs a pseudo version of the secret, however, the copy-pasted secret is stored_) 
# MAGIC ```
# MAGIC
# MAGIC databricks secrets get-secret {your_scope_name_prefix}-azure-openai openai_api_key 
# MAGIC  
# MAGIC ```
# MAGIC
# MAGIC
# MAGIC Likewise, we can check the defined scope and secret using `dbuitls.secrets.get(scope, key)` -- it will be shown as `['REDACTED]'` 
# MAGIC
# MAGIC ```
# MAGIC  
# MAGIC dbutils.secrets.get("{your_scope_name_prefix}-azure-openai", "openai_api_key")
# MAGIC
# MAGIC '[REDACTED]'
# MAGIC ```
# MAGIC
# MAGIC

# COMMAND ----------

# DBTITLE 1,Create Scope with Databricks SDK
## If you don't already have one created ...

from databricks.sdk import WorkspaceClient

w = WorkspaceClient()

scope_name= "<{your_scope_name_prefix}_azure_openai>"  ## update to use <ServicePrinciple> reference
scope_key= "<azopenai_api_key>"

w.secrets.create_scope(scope_name)

# COMMAND ----------

# DBTITLE 1,Create Scope Key with Secret
w.secrets.put_secret(scope_name,scope_key,string_value ="<secret>") ## do not leave secret exposed 

# COMMAND ----------

# DBTITLE 1,Read/Get_secret
# w.secrets.get_secret(scope_name,scope_key).value
# (NB: this outputs a pseudo version of the secret, however, the copy-pasted secret is stored)

# COMMAND ----------

# MAGIC %md
# MAGIC ### [3] Define & Serve the `Azure Openai GPT-4o` External Model as an Endpoint 

# COMMAND ----------

# DBTITLE 1,run notebook utils
# MAGIC %run ./utils

# COMMAND ----------

# DBTITLE 1,UC variables
remove_widgets() 
uc_config = setup_uc_paths(spark=None, use_widgets=False); ## if you update the values in widgets -- it will automatically trigger an update of the UC paths

# Extract catalog, schema, volume names
catalog_name = uc_config["catalog_name"]
schema_name = uc_config["schema_name"]
volume_name = uc_config["volume_name"]
# external_endpoint_name = uc_config["external_endpoint_name"]


# COMMAND ----------

# DBTITLE 1,Set the environment variables
import os

os.environ["DATABRICKS_HOST"] = "https://e2-demo-west.cloud.databricks.com/"
os.environ["DATABRICKS_TOKEN"] = dbutils.secrets.get("mmt", "databricks_token") ## PAT token | Best Practice: Service Principal (SP) token is recommended for production; require privileges/permissions to create SP.

## Set the environment variables
# os.environ["DATABRICKS_HOST"] = "https://{workspace-instance}.{cloud or shard}.databricks.com/"
# os.environ["DATABRICKS_TOKEN"] = dbutils.secrets.get({scope_name}, {scope_key}) ## PAT token | Best Practice: Service Principal (SP) token is recommended for production; require privileges/permissions to create SP. 

# COMMAND ----------

# DBTITLE 1,Use serving-endpoints API
import os
import requests
import json

# Retrieve the Databricks host and token from environment variables
databricks_host = os.getenv("DATABRICKS_HOST")
databricks_token = os.getenv("DATABRICKS_TOKEN")

# Check if the environment variables are set
if not databricks_host or not databricks_token:
    raise ValueError("Databricks host and token must be set in environment variables")

# Define the endpoint URL
url = f"{databricks_host}/api/2.0/serving-endpoints"

# Define the headers
headers = {
    "Authorization": f"Bearer {databricks_token}",
    "Content-Type": "application/json"
}

# Define the payload
payload = {
    "name": "<external_endpoint_name e.g. az_openai_gpt4o>",
    "config": {
        "served_entities": [
            {
                "name": "az-openai-completions",                            

                "external_model": {
                    "name": "gpt-4o",
                    "provider": "openai",
                    "task": "llm/v1/chat",
                    "openai_config": {
                        "openai_api_type": "azure",                        
                        "openai_api_key": f"{{{{secrets/{scope_name}/{scope_key}}}}}",                        
                        "openai_api_base": "https://<served_ai_foundry_endpoint_name>.openai.azure.com/",                         
                        # "openai_api_base": "https://hls-fedemo-azure-openai.openai.azure.com/",
                        # "openai_api_key": "{{secrets/hls_fedemo_azure_openai/azopenai_api_key}}", #databricks cli/sdk registered 
                        
                        "openai_deployment_name": "gpt-4o-2024-11-20", 
                        "openai_api_version": "2025-01-01-preview"
                    }
                }
            }
        ],
    },        
    "ai_gateway": {
      "usage_tracking_config": {
        "enabled": True
      },
      "inference_table_config": {
          "catalog_name": "demos_genie",
          "schema_name": "hls_ai_drug_discovery",
          "table_name_prefix": "", #"your_table_prefix",
          "enabled": True
      }
    },
      "tags": [
          {
              "key": "removeAfter",
              "value": "2026-01-31"
          },
          {
              "key": "<project>",
              "value": "<name_of_project e.g. ai-driven_drug_discovery>"
          },
          {
              "key": "do-not-delete",
              "value": "True"
          }
      ]
  }

# Make the POST request
response = requests.post(url, headers=headers, data=json.dumps(payload))

# Check the response
if response.status_code == 200:
    print("Endpoint created successfully")
else:
    print(f"Failed to create endpoint: {response.status_code}")
    print(response.text)

# COMMAND ----------

# DBTITLE 1,Deployed Serving endpoint
# MAGIC %md
# MAGIC #### Deployed serving endpoint Info.:
# MAGIC - ws endpoint: [az_openai_gpt4o](https://e2-demo-west.cloud.databricks.com/ml/endpoints/az_openai_gpt4o?o=2556758628403379)
# MAGIC - serving-endpoint: `https://e2-demo-west.cloud.databricks.com/serving-endpoints/az_openai_gpt4o/invocations`
# MAGIC - [az_openai_gpt4o_payload](https://e2-demo-west.cloud.databricks.com/explore/data/demos_genie/hls_ai_drug_discovery/az_openai_gpt4o_payload?o=2556758628403379): `demos_genie.hls_ai_drug_discovery.az_openai_gpt4o_payload`

# COMMAND ----------

# MAGIC %md
# MAGIC ### [4] Test using the External Endpoint for inferencing via [`ai_query()`](https://docs.databricks.com/en/large-language-models/ai-functions.html#ai_query)
# MAGIC
# MAGIC Now that the endpoint has been successfully created, we can test making inferencing before using it in our AI/BI & Genie Space.       
# MAGIC Here we will demonstrate batch inferencing using `ai_query()`.    
# MAGIC
# MAGIC Below is a template of how we setup this batch inference in SQL:
# MAGIC
# MAGIC ```
# MAGIC
# MAGIC SELECT
# MAGIC   {input_column},   -- Placeholder for the input column
# MAGIC   ai_query(
# MAGIC     'az_openai_gpt4o',
# MAGIC     CONCAT({prompt_placeholder}, {input_column})    -- Placeholder for the prompt and input
# MAGIC   ) AS {output_column}  -- Placeholder for the output column
# MAGIC FROM {table_name};  -- Placeholder for the table name
# MAGIC LIMIT 50
# MAGIC ```
# MAGIC
# MAGIC - Note: [Pay-per-token foundation models](https://docs.databricks.com/en/machine-learning/foundation-model-apis/index.html#pay-per-token-foundation-model-apis) are good for testing but have limits.    
# MAGIC - They are unsuitable for processing more than 100 rows; use [provisioned foundation models](https://docs.databricks.com/en/machine-learning/foundation-model-apis/deploy-prov-throughput-foundation-model-apis.html) for larger datasets.
# MAGIC

# COMMAND ----------

# DBTITLE 1,TEST external gpt4o endpoint
from pyspark.sql import functions as F

# Parameters
# catalog_name = "mmt_demos2" #"<your_catalog_name>" 
# schema_name = "ai_driven_drug_discovery" #"<your_schema_name>"
# endpoint_name = "az_openai_gpt4o" ## name of deployed external endpoint on your workspace
score_threshold = 0.85
organism_filter = "%human%"

# Read tables
protein_df = spark.table(f"{catalog_name}.{schema_name}.proteinclassification_tiny")
organism_df = spark.table(f"{catalog_name}.{schema_name}.tinysample_organism_info_scientificNsimple")

# Filter and join
filtered_df = (
    protein_df
    .filter(F.col("score") > score_threshold)
    .join(
        organism_df,
        protein_df.OrganismName == organism_df.OrganismName,
        "inner"
    )
    .filter(F.lower(F.col("Organism_SimpleTerm")).like(organism_filter.lower()))
)

# Create AI query prompt
prompt_template = (
    "You are well-versed in membrane proteins and drug discovery research. "
    "Please be brief. Provide a dictionary of responses in the format "
    '{"key": "response"} to the following keys "information", "recent_research" '
    'and highlight "under_researched_areas" that hold promise for drug discovery '
    'for the given "protein name": {protein_name} '
    'Output just the dictionary {"information": "response", "recent_research": "response", '
    '"under_researched_areas": "response"}. Do not include "```json" strings in output'
)

# Build final dataframe with AI query
result_df = (
    filtered_df
    .withColumn(
        "ai_prompt",
        F.concat(
            F.lit(prompt_template.replace("{protein_name}", "")),
            F.col("ProteinName")
        )
    )
    .withColumn(
        "researchDict",
        F.expr(f"ai_query('{endpoint_name}', ai_prompt)")
    )
    .withColumn(
        "information",
        F.get_json_object(F.col("researchDict"), "$.information")
    )
    .withColumn(
        "recent_research",
        F.get_json_object(F.col("researchDict"), "$.recent_research")
    )
    .withColumn(
        "under_researched_areas",
        F.get_json_object(F.col("researchDict"), "$.under_researched_areas")
    )
    .select(
        organism_df.OrganismName,
        F.col("Organism_SimpleTerm"),
        protein_df.ProteinName,
        F.col("researchDict"),
        F.col("information"),
        F.col("recent_research"),
        F.col("under_researched_areas"),
        F.col("label").alias("ProteinType"),
        F.col("score").alias("ProteinClassificationScore")
    )
)

# Display or save results
display(result_df)

# COMMAND ----------

# MAGIC %md
# MAGIC ### [5] Register the [`ai_query()`](https://docs.databricks.com/aws/en/sql/language-manual/functions/ai_query) as `SQL_function` to Unity Catalog 
# MAGIC Similarly to previous example, We can register the `ai_query()` calling the external Foundation Model `az_openai_gpt4o` to help with getting protein related research infomation associated with Organism of interest. This makes it easier for calling the registered SQL_function later in either AI/BI Dashboard/Genie Space or even in the Playground. 

# COMMAND ----------

# DBTITLE 1,[x] Register as ai_query() as SQL function
### via PySpark 
# Define the SQL function

create_function_query = f"""
CREATE OR REPLACE FUNCTION {catalog_name}.{schema_name}.get_protein_research_info(ProteinName STRING)
RETURNS STRING
COMMENT 'Returns a dictionary of information, recent_research, and under_researched_areas as JSON output.'
RETURN ai_query(
          'az_openai_gpt4o',
          CONCAT(
            'You are well-versed in membrane proteins and drug discovery research. Please be brief. Provide a dictionary of responses in the format {{"key": "response"}} to the following keys "information", "recent_research" and highlight "under_researched_areas" that hold promise for drug discovery for the given "protein name": ',
            ProteinName,
            '. Output just the dictionary {{"information": "response", "recent_research": "response", "under_researched_areas": "response"}}. Do not include ```json" strings in output.')
);
"""

# Execute the query to create the function
spark.sql(create_function_query)

# COMMAND ----------

# DBTITLE 1,Register as ai_query() as SQL function
# Configuration
# AI_MODEL = "<external_endpoint_name e.g. az_openai_gpt4o>"
AI_MODEL = external_endpoint_name
FUNCTION_NAME = f"{catalog_name}.{schema_name}.get_protein_research_info"

# System prompt (easier to maintain)
SYSTEM_PROMPT = """You are a membrane proteins and drug discovery expert.
Analyze the given protein and return ONLY valid JSON with these exact keys:
- information: brief protein overview
- recent_research: key recent findings  
- under_researched_areas: promising drug discovery opportunities

Output format: {"information": "...", "recent_research": "...", "under_researched_areas": "..."}
No markdown. No code blocks. Be concise."""

create_function_query = f"""
CREATE OR REPLACE FUNCTION {FUNCTION_NAME}(ProteinName STRING)
RETURNS STRUCT<information: STRING, recent_research: STRING, under_researched_areas: STRING>
RETURN 
  FROM_JSON(
    ai_query(
      '{AI_MODEL}',
      CONCAT('{SYSTEM_PROMPT}\n\nProtein: "', REPLACE(ProteinName, '"', ''), '"')
    ),
    'STRUCT<information: STRING, recent_research: STRING, under_researched_areas: STRING>'
  );
"""

spark.sql(create_function_query)

# COMMAND ----------

# DBTITLE 1,Test Bulk process with SQL_function
# Parameters

score_threshold = 0.8
organism_filter = "%Zebrafish%"  # Options: "%human%", "%Zebrafish%", etc.
output_table_name = "tinysample_organism_protein_research_info"
save_results = False  # Set to True to save results

# Build the query
query = f"""
SELECT 
  OrganismName,
  Organism_SimpleTerm,
  ProteinName,
  researchDict.information AS information,
  researchDict.recent_research AS recent_research,
  researchDict.under_researched_areas AS under_researched_areas,
  ProteinType,
  ProteinClassificationScore
FROM (
  SELECT
    p.ProteinName,
    p.label AS ProteinType,
    p.score AS ProteinClassificationScore,
    p.GeneName,
    p.Sequence,
    p.Molecular_Weight,
    p.OrganismName,
    o.Organism_SimpleTerm,
    {catalog_name}.{schema_name}.get_protein_research_info(p.ProteinName) AS researchDict
  FROM
    {catalog_name}.{schema_name}.proteinclassification_tiny p
  JOIN 
    {catalog_name}.{schema_name}.tinysample_organism_info_scientificNsimple o 
    ON p.OrganismName = o.OrganismName
  WHERE
    p.score > {score_threshold}
    AND o.Organism_SimpleTerm ILIKE '{organism_filter}'
)
"""

# Display the query
print("Generated SQL Query:")
print("=" * 80)
print(query)
print("=" * 80)
print(f"\nQuery Parameters:")
print(f"  - Catalog: {catalog_name}")
print(f"  - Schema: {schema_name}")
print(f"  - Score Threshold: {score_threshold}")
print(f"  - Organism Filter: {organism_filter}")
print("=" * 80)

# Execute the query
print("\nExecuting query...")
sDF_protein_research_info = spark.sql(query)

# Display results
print(f"\nQuery executed successfully. Displaying results...")
display(sDF_protein_research_info)

# Optional: Save results to table
if save_results:
    output_table = f"{catalog_name}.{schema_name}.{output_table_name}"
    print(f"\nSaving results to table: {output_table}")
    sDF_protein_research_info.write.mode("overwrite").option("mergeSchema", "true").saveAsTable(output_table)
    print(f"Results saved successfully to {output_table}")
else:
    print("\nNote: Results not saved. Set save_results=True to save to table.")

# COMMAND ----------

# MAGIC %md
# MAGIC ### [6] Bulk process to derive `protein_research_info` 
# MAGIC - In order to more efficiently explore without waiting for the queries to run for specific `Organism_SimpleTerm` and `ProteinClassificationScore` during ad-hoc explorations, we will pre-process to extract the information pertaining to `recent research` and `"under_researched_areas" that hold promise for drug discovery` for the list of proteins in our `tinysample`. 
# MAGIC - Although this will take a bit of time up-front, it will allow for more efficient and simpler filtering of the inferred proteins and the research-relevant outputs. 
# MAGIC - For production use, [provisioned throughput model serving endpoint](https://docs.databricks.com/en/machine-learning/foundation-model-apis/deploy-prov-throughput-foundation-model-apis.html#provisioned-throughput-endpoint-ui) is highly recommended if [Foundation Model APIs](https://docs.databricks.com/aws/en/machine-learning/model-serving/model-serving-limits#provisioned-throughput-limits) can be leveraged!

# COMMAND ----------

# DBTITLE 1,[x] pyspark Bulk process with SQL_function
# from pyspark.sql.functions import col, get_json_object, expr

# # Define the base query
# base_query = f"""
# SELECT
#   p.ProteinName,
#   p.label AS ProteinType,
#   p.score AS ProteinClassificationScore,
#   p.GeneName,
#   p.Sequence,
#   p.Molecular_Weight,
#   p.OrganismName,
#   o.Organism_SimpleTerm
# FROM
#    {catalog_name}.{schema_name}.proteinclassification_tiny p
#   JOIN  {catalog_name}.{schema_name}.tinysample_organism_info_scientificNsimple o 
#   ON p.OrganismName = o.OrganismName
# """

# # Create a DataFrame from the base query
# df_base = spark.sql(base_query)

# # Apply the get_protein_research_info function and extract JSON fields
# df_research_info = df_base.withColumn(
#     "researchDict", 
#     expr(f"{catalog_name}.{schema_name}.get_protein_research_info(ProteinName)").cast('string')
# ).withColumn(
#     "information", 
#     get_json_object(col("researchDict"), "$.information")
# ).withColumn(
#     "recent_research", 
#     get_json_object(col("researchDict"), "$.recent_research")
# ).withColumn(
#     "under_researched_areas", 
#     get_json_object(col("researchDict"), "$.under_researched_areas")
# ).select(
#     "OrganismName",
#     "Organism_SimpleTerm",
#     "ProteinName",
#     "researchDict",
#     "information",
#     "recent_research",
#     "under_researched_areas",
#     "ProteinType",
#     "ProteinClassificationScore"
# )


# # df_base.count(), df_research_info.count()

# # # Write the DataFrame to a Delta Table 
# # df_research_info.write.mode("overwrite").option("mergeSchema", "true").saveAsTable(f"{catalog_name}.{schema_name}.tinysample_organism_protein_research_info")

# COMMAND ----------

# DBTITLE 1,pyspark Bulk process with SQL_function
from pyspark.sql import functions as F

# Parameters

score_threshold = None  # Options: None (no filter), 0.8, 0.9, etc.
organism_filter = None  # Options: None (all organisms), "%human%", "%Zebrafish%", etc.
output_table_name = "tinysample_organism_protein_research_info"
save_results = True  # Set to False to skip saving
display_results = True  # Set to False to skip displaying results
display_limit = 100  # Options: None (show all), 100, 500, etc.

# Build the WHERE clause conditionally
where_conditions = []

if score_threshold is not None:
    where_conditions.append(f"p.score > {score_threshold}")

if organism_filter is not None:
    where_conditions.append(f"o.Organism_SimpleTerm ILIKE '{organism_filter}'")

# Only add WHERE clause if there are conditions
where_clause = f"WHERE {' AND '.join(where_conditions)}" if where_conditions else ""

# Build the base query
base_query = f"""
SELECT
  p.ProteinName,
  p.label AS ProteinType,
  p.score AS ProteinClassificationScore,
  p.GeneName,
  p.Sequence,
  p.Molecular_Weight,
  p.OrganismName,
  o.Organism_SimpleTerm
FROM
  {catalog_name}.{schema_name}.proteinclassification_tiny p
JOIN 
  {catalog_name}.{schema_name}.tinysample_organism_info_scientificNsimple o 
  ON p.OrganismName = o.OrganismName
{where_clause}
"""

# Display the query
print("Generated Base SQL Query:")
print("=" * 80)
print(base_query)
print("=" * 80)
print(f"\nQuery Parameters:")
print(f"  - Catalog: {catalog_name}")
print(f"  - Schema: {schema_name}")
print(f"  - Score Threshold: {score_threshold if score_threshold is not None else 'None (no filter)'}")
print(f"  - Organism Filter: {organism_filter if organism_filter else 'None (no filter)'}")
print(f"  - Display Results: {display_results}")
print(f"  - Display Limit: {display_limit if display_limit else 'None (show all)'}")
print("=" * 80)

# Execute base query
print("\nExecuting base query...")
df_base = spark.sql(base_query)
print(f"Base query executed. Row count: {df_base.count()}")

# Apply the get_protein_research_info function and extract STRUCT fields
print("\nApplying protein research info function and extracting STRUCT fields...")
df_research_info = (
    df_base
    .withColumn(
        "researchDict", 
        F.expr(f"{catalog_name}.{schema_name}.get_protein_research_info(ProteinName)")
    )
    .withColumn(
        "information", 
        F.col("researchDict.information")
    )
    .withColumn(
        "recent_research", 
        F.col("researchDict.recent_research")
    )
    .withColumn(
        "under_researched_areas", 
        F.col("researchDict.under_researched_areas")
    )
    .select(
        "OrganismName",
        "Organism_SimpleTerm",
        "ProteinName",
        "researchDict",
        "information",
        "recent_research",
        "under_researched_areas",
        "ProteinType",
        "ProteinClassificationScore"
    )
)

# Display results (conditionally)
if display_results:
    print("\nDisplaying research info results...")
    if display_limit is not None:
        print(f"(Limited to {display_limit} rows)")
        display(df_research_info.limit(display_limit))
    else:
        print("(Showing all rows)")
        display(df_research_info)
else:
    print("\nNote: Display skipped. Set display_results=True to display results.")

# Save to Delta Table
if save_results:
    output_table = f"{catalog_name}.{schema_name}.{output_table_name}"
    print(f"\nSaving results to table: {output_table}")
    
    # Drop table first to avoid schema conflicts
    spark.sql(f"DROP TABLE IF EXISTS {output_table}")
    
    # Save without mergeSchema option
    df_research_info.write.mode("overwrite").saveAsTable(output_table)
    print(f"Results saved successfully to {output_table}")
else:
    print("\nNote: Results not saved. Set save_results=True to save to table.")

# Summary
print("\n" + "=" * 80)
print("EXECUTION SUMMARY")
print("=" * 80)
print(f"Total proteins processed: {df_research_info.count()}")
print(f"Score threshold: {score_threshold if score_threshold is not None else 'None (no filter)'}")
print(f"Organism filter applied: {organism_filter if organism_filter else 'None (no filter)'}")
print(f"Results displayed: {display_results}")
if display_results and display_limit:
    print(f"Display limit: {display_limit} rows")
print(f"Results saved: {save_results}")
if save_results:
    print(f"Output table: {catalog_name}.{schema_name}.{output_table_name}")
print("=" * 80)

# COMMAND ----------

# DBTITLE 1,Check/Read tinysample_organism_protein_research_info
sDF_protein_research_info = spark.table(f"{catalog_name}.{schema_name}.tinysample_organism_protein_research_info")  

display(sDF_protein_research_info)

# COMMAND ----------


