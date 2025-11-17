# AI-Driven Drug Discovery Notebooks:

This set of notebooks found within this repository folder involves data preparation and foundational model endpoint setup and usage.  

Most of the notebooks will rely on the [**`utils.py`**](notebooks/utils.py) to configure the Unity Catalog paths.    
**Users are required to update the Unity Catalog `<catalog_name>`** :

```
# ===============================================================
# REQUIRED: UPDATE THESE PLACEHOLDER VALUES FOR YOUR ENVIRONMENT
# ===============================================================
# These values are used as placeholder / defaults 

## REPLACE THE VALUES BELOW WITH YOUR OWN BEFORE running the notebooks:

CATALOG_NAME = "<your_catalog_name>"      # TODO: Replace with <your_catalog_name>
SCHEMA_NAME = "ai_driven_drug_discovery"  # TODO: Use Default OR Replace with <your_schema_name>
VOLUME_NAME = "protein_seq"               # TODO: Use Default OR Replace with <your_volume_name> for file storage
ENDPOINT_NAME = "az_openai_gpt4o"         # TODO: Use Default OR Replace with <your_external_endpoint_name> e.g. AI Gateway endpoint name
    
```    

### The notebooks step through the workflow: 

1. Programatically [`Download_UNIPROT_fasta`](notebooks/1.0_Download_UNIPROT_fasta.ipynb) to Unity Catalog Volumes 

2. Run the [`ProteinData_ETL`](notebooks/2.0_ProteinData_ETL.ipynb) to clean and enrich the raw fasta data e.g. with molecular weights info.
    - For this ETL pipeline notebook, users are required to update `<catalog_name>` within before using it as a source for setting up, validating, and running the notebook as an ETL pipeline. 

3. Apply [`TransformerCNN_Protein_Classification`](notebooks/3.0_TransformerCNN_Protein_Classification.ipynb) to assess if protein sequence is likely of the type `water-soluble` or `cell-membrane transport` function specific. 

4. Leverage Databricks's Foundational LLMs (FMAPI) and AI Functions to [`Register&Query_scientific2simple_UDF`](notebooks/4.0_Register&Query_scientific2simple_UDF.ipynb) e.g. useful for democratizing scientific knowledge and can be applied in bulk data processing using `ai_query()`.  

5. Serve and harness [`ExternalEndpoint`](notebooks/5.0_Create&Query_ExternalEndpoint_gpt4o.ipynb) on Databricks Model Serving endpoint to facilitate downstream discovery efforts.    

<br> 

### Downstream: Setup up Dashboard

Once these notebooks are run, we can configure the AIBI Dashboard with Genie using the associate template within [`./dashboards`](../dashboards)  folder.

<br>    

---   