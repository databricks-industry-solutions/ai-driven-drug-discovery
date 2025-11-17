# Databricks notebook source
# DBTITLE 1,Please specify UC config before starting
# utils.py | .ipynb
# ===============================================================
# REQUIRED: UPDATE THESE PLACEHOLDER VALUES FOR YOUR ENVIRONMENT
# ===============================================================
# These values are used as placeholder / defaults 

## REPLACE THE VALUES BELOW WITH YOUR OWN BEFORE running the notebooks:

CATALOG_NAME = "<your_catalog_name>"      # TODO: Replace with <your_catalog_name>
SCHEMA_NAME = "ai_driven_drug_discovery"  # TODO: Use Default OR Replace with <your_schema_name>
VOLUME_NAME = "protein_seq"               # TODO: Use Default OR Replace with <your_volume_name> for file storage
ENDPOINT_NAME = "az_openai_gpt4o"         # TODO: Use Default OR Replace with <your_external_endpoint_name> e.g. AI Gateway endpoint name

from pyspark.sql import SparkSession

# COMMAND ----------

# DBTITLE 1,set up  UC paths configs
def setup_uc_paths(spark: SparkSession = None, 
                     use_widgets: bool = True) -> dict:
    """
    Setup Unity Catalog resources and return configuration
    
    Args:
        spark: SparkSession (optional - auto-detected if None)
        use_widgets: If True (default), creates widgets for configuration
    
    Returns:
        Dictionary with catalog_name, schema_name, volume_name, external_endpoint_name, 
        volume_location, schema_path, volume_path
    
    Example:
        config = setup_environment()
        print(config['catalog_name'])
        print(config['volume_location'])
    """
    # Auto-detect spark session if not provided
    if spark is None:
        spark = SparkSession.getActiveSession()
        if spark is None:
            raise RuntimeError(
                "No active Spark session found. "
                "This should not happen in Databricks notebooks."
            )
    
    # Setup widgets if requested
    if use_widgets:
        try:
            dbutils.widgets.text("catalog_name", CATALOG_NAME, "1. Catalog Name")
            dbutils.widgets.text("schema_name", SCHEMA_NAME, "2. Schema Name")
            dbutils.widgets.text("volume_name", VOLUME_NAME, "3. Volume Name")
            dbutils.widgets.text("external_endpoint_name", ENDPOINT_NAME, "4. AI Gateway Endpoint")
        except:
            pass  # Widgets may already exist
        
        # Read values from widgets, fallback to defaults
        catalog_name = dbutils.widgets.get("catalog_name") or CATALOG_NAME
        schema_name = dbutils.widgets.get("schema_name") or SCHEMA_NAME
        volume_name = dbutils.widgets.get("volume_name") or VOLUME_NAME
        external_endpoint_name = dbutils.widgets.get("external_endpoint_name") or ENDPOINT_NAME
    else:
        # Use only hardcoded defaults
        catalog_name = CATALOG_NAME
        schema_name = SCHEMA_NAME
        volume_name = VOLUME_NAME
        external_endpoint_name = ENDPOINT_NAME
    
    # Check for placeholder values
    placeholder_configs = []
    if catalog_name.startswith("<") and catalog_name.endswith(">"):
        placeholder_configs.append("Catalog Name")
    if schema_name.startswith("<") and schema_name.endswith(">"):
        placeholder_configs.append("Schema Name")
    if volume_name.startswith("<") and volume_name.endswith(">"):
        placeholder_configs.append("Volume Name")
    if external_endpoint_name.startswith("<") and external_endpoint_name.endswith(">"):
        placeholder_configs.append("External Endpoint Name")
    
    if placeholder_configs:
        if use_widgets:
            raise ValueError(
                f"Placeholder values detected. Please fill in the following widget(s) at the top of the notebook:\n"
                f"  - {', '.join(placeholder_configs)}\n"
                f"Then re-run this cell."
            )
        else:
            raise ValueError(
                f"Placeholder values detected. Please update the following value(s) at the top of utils.py:\n"
                f"  - {', '.join(placeholder_configs)}\n"
                f"Replace the placeholder values (e.g., '<your_catalog_name>') with actual values."
            )
    
    # Validate required configuration
    missing_configs = []
    if not catalog_name:
        missing_configs.append("Catalog Name")
    if not schema_name:
        missing_configs.append("Schema Name")
    if not volume_name:
        missing_configs.append("Volume Name")
    if not external_endpoint_name:
        missing_configs.append("External Endpoint Name")
    
    if missing_configs:
        if use_widgets:
            raise ValueError(
                f"Please fill in the following widget(s) at the top of the notebook:\n"
                f"  - {', '.join(missing_configs)}\n"
                f"Then re-run this cell."
            )
        else:
            raise ValueError(
                f"Please set the following value(s) at the top of utils.py:\n"
                f"  - {', '.join(missing_configs)}"
            )
    
    # Create UC resources if they don't exist
    spark.sql(f"CREATE CATALOG IF NOT EXISTS {catalog_name}")
    spark.sql(f"CREATE SCHEMA IF NOT EXISTS {catalog_name}.{schema_name}")
    spark.sql(f"CREATE VOLUME IF NOT EXISTS {catalog_name}.{schema_name}.{volume_name}")
    
    # Calculate derived paths
    volume_location = f"/Volumes/{catalog_name}/{schema_name}/{volume_name}"
    schema_path = f"{catalog_name}.{schema_name}"
    volume_path = f"{catalog_name}.{schema_name}.{volume_name}"

    uc_config = {
        'catalog_name': catalog_name,
        'schema_name': schema_name,
        'volume_name': volume_name,
        'external_endpoint_name': external_endpoint_name,
        'volume_location': volume_location,
        'schema_path': schema_path,
        'volume_path': volume_path
    } 
    
    # Print summary
    config_source = " (from widgets)" if use_widgets else " (from defaults)"
    print("="*70)
    print(f"UC Paths Configured{config_source}")
    print("="*70)
    for key, value in uc_config.items():
        print(f"{key}: {value}")
    print("="*70)
    
    return uc_config


def remove_widgets() -> None:
    """
    Remove configuration widgets from notebook UI
    """
    try:
        dbutils.widgets.removeAll()
        print("Existing Widgets Removed")
    except:
        print("No widgets to remove")

# COMMAND ----------

# DBTITLE 1,test/ example usage
# remove_widgets() 
# setup_uc_paths(spark=None, use_widgets=False);

# COMMAND ----------

# MAGIC %md
# MAGIC
