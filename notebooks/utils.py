# Databricks notebook source
# DBTITLE 1,Please specify UC config before starting
# utils.py | .ipynb
# ===============================================================
# REQUIRED: UPDATE THESE PLACEHOLDER VALUES FOR YOUR ENVIRONMENT
# ===============================================================
# These values are used as placeholder / defaults 

## REPLACE THE VALUES BELOW WITH YOUR OWN BEFORE running the notebooks:

CATALOG_NAME = "mmt_demos2"
# CATALOG_NAME = "<your_catalog_name>"      # TODO: Replace with <your_catalog_name>
SCHEMA_NAME = "ai_driven_drug_discovery"  # TODO: Use Default OR Replace with <your_schema_name>
VOLUME_NAME = "protein_seq"               # TODO: Use Default OR Replace with <your_volume_name> for file storage
ENDPOINT_NAME = "az_openai_gpt4o"         # TODO: Use Default OR Replace with <your_external_endpoint_name> e.g. AI Gateway endpoint name

# COMMAND ----------

# DBTITLE 1,set up  UC paths configs
def setup_uc_paths(use_widgets: bool = True, print_endpoint: bool = False, silent: bool = True) -> dict:
    """
    Setup Unity Catalog resources and return configuration
    
    Args:
        use_widgets: If True (default), creates widgets for configuration
        print_endpoint: If True, includes external_endpoint_name in printed output (default: False)
        silent: If True (default), suppresses all printed output
    
    Returns:
        Dictionary with catalog_name, schema_name, volume_name, external_endpoint_name, 
        volume_location, schema_path, volume_path
    
    Example:
        # Silent mode (default - no output)
        config = setup_uc_paths()
        
        # With output printed (no endpoint)
        config = setup_uc_paths(silent=False)
        
        # With output and endpoint printed
        config = setup_uc_paths(silent=False, print_endpoint=True)
    """
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
    
    # Create UC resources if they don't exist (using global spark object)
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
    
    # Print summary with selective fields (unless silent mode)
    if not silent:
        config_source = " (from widgets)" if use_widgets else " (from defaults)"
        print("="*70)
        print(f"UC Paths Configured{config_source}")
        print("="*70)
        for key, value in uc_config.items():
            # Skip external_endpoint_name if print_endpoint is False
            if key == 'external_endpoint_name' and not print_endpoint:
                continue
            print(f"{key}: {value}")
        print("="*70)
    
    return uc_config


def remove_widgets(silent: bool = True) -> None:
    """
    Remove configuration widgets from notebook UI
    
    Args:
        silent: If True (default), suppresses printed output
    """
    try:
        dbutils.widgets.removeAll()
        if not silent:
            print("Existing Widgets Removed")
    except:
        if not silent:
            print("No widgets to remove")

# COMMAND ----------

# DBTITLE 1,test/ example usage
# setup_uc_paths(use_widgets=False, print_endpoint=False, silent=True)

# Test: With output printed (but no endpoint)
# setup_uc_paths(use_widgets=False, print_endpoint=False, silent=False)

# COMMAND ----------

# MAGIC %md
# MAGIC
