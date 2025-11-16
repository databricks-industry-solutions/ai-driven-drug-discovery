# Databricks notebook source
# MAGIC %md
# MAGIC # Protein Classification 
# MAGIC ## Water `Soluble` vs `Membrane` Transport Proteins

# COMMAND ----------

# MAGIC %md
# MAGIC ### Proteins: Form & Function
# MAGIC
# MAGIC Proteins are often considered the "_working molecules_" within every living organism as they perform a myriad of vital functions, including providing structural support, catalyzing chemical reactions as enzymes, transporting molecules, acting as hormones, defending against infection, and playing a key role in cell signaling and division ([ref](https://en.wikipedia.org/wiki/Protein)). They can exist in different shapes, sizes, alone or as part of a multi-unit structure, and they can alter their shape frequently or remain virtually static. Their structural diversity arise from their unique (genetically determined) composition of amino acid sequence. When fully folded, their distinct __surface characteristics__ determine which other molecules they interact with, e.g. their conformation can change in subtle or dramatic ways when they bind with other molecules ([ref](https://www.nature.com/scitable/topicpage/protein-function-14123348/)).    
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC <!-- [maybe insert image here]   -->
# MAGIC <!-- <img src="https://storage.googleapis.com/gdm-deepmind-com-prod-public/media/original_images/62273f8719ed3b2a84c8fd13_Fig201.svg"/> -->
# MAGIC <img src="https://storage.googleapis.com/gdm-deepmind-com-prod-public/media/original_images/62273f8719ed3b2a84c8fd13_Fig201.svg" style="padding: 10px;" alt="Protein Structure" title="Complex 3D shapes emerge from a string of amino acids. Source:https://deepmind.google/discover/blog/alphafold-using-ai-for-scientific-discovery-2020/" width="1000"/> 
# MAGIC
# MAGIC
# MAGIC `Source:`[`https://deepmind.google/discover/blog/alphafold-using-ai-for-scientific-discovery-2020`](https://deepmind.google/discover/blog/alphafold-using-ai-for-scientific-discovery-2020)
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC Proteins can generally be categorized as either water-`soluble` proteins or `membrane` transport proteins, based on their amino acid composition: `Soluble` proteins have a hydrophilic surface, allowing them to interact with water molecules and dissolve readily in the cell cytoplasm, while `membrane` transport proteins have hydrophobic regions that interact with the lipid bilayer of the cell membrane, allowing them to embed within the membrane ([ref](https://www.nature.com/articles/s41598-017-17216-1)).    

# COMMAND ----------

# MAGIC %md
# MAGIC ##### Example of membrane transport proteins     
# MAGIC ###          
# MAGIC
# MAGIC <img src= "https://ars.els-cdn.com/content/image/3-s2.0-B9780128042540000016-f01-03-9780128042540.jpg" width="600">
# MAGIC
# MAGIC `Source:`[`Carlson, B.M. in The Human Body, 2019, Chapter 1, Cells`](https://www.sciencedirect.com/topics/neuroscience/membrane-channel)
# MAGIC
# MAGIC
# MAGIC ---    
# MAGIC ##### Example of how membrane transport protein is visualized    
# MAGIC
# MAGIC
# MAGIC <img src="https://pub.mdpi-res.com/cancers/cancers-12-01624/article_deploy/html/images/cancers-12-01624-g002.png?1592552138" width="800">
# MAGIC
# MAGIC ```
# MAGIC Na+-, K+-ATPase overall structure.    
# MAGIC
# MAGIC (A) Na+-, K+-ATPase consists of a catalytic α subunit (alpha helices) and a regulatory β subunit (beta pleated sheets). 
# MAGIC The α subunit consists of 10 transmembrane helices, harboring 3 different cytoplasmic domains: 
# MAGIC the actuator responsible for dephosphorylation (shown in red); the nucleotide-binding, responsible for ATP binding (shown in blue); 
# MAGIC and the phosphorylation domains (shown in cyan). 
# MAGIC The β subunit consists of one transmembrane helix with a large glycosylated extracellular domain (shown in hexagon boxes). 
# MAGIC ECM = extracellular milieu; CYT = cytoplasm.    
# MAGIC
# MAGIC (B) Overall domain architecture of Na+/K+ transporter in the Na+-bound state (Protein Data Bank [PDB] code 4HQJ). 
# MAGIC Catalytic α subunit is colored in blue, β subunit is shown in yellow, and Na+ ions are shown in red.
# MAGIC ```         
# MAGIC `Source:` [```Almasi S, El Hiani Y. Exploring the Therapeutic Potential of Membrane Transport Proteins:    
# MAGIC Focus on Cancer and Chemoresistance. Cancers. 2020; 12(6):1624. https://doi.org/10.3390/cancers12061624```](https://doi.org/10.3390/cancers12061624)
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC **The ability to classify protein sequences into functional types** ([ref](https://portlandpress.com/essaysbiochem/article/66/3/255/231650/Uncovering-protein-function-from-classification-to)) e.g. such as _water-soluble_ or _membrane transport type_, **can help provide crucial information about their potential location and function within a cell**, allowing researchers to better understand how a protein might interact with other molecules and contribute to cellular processes, particularly regarding the transport of substances across the cell membrane. (refs: e.g. [i](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3275-6),[ii](https://www.nature.com/articles/s41598-017-17216-1),[iii](https://www.nature.com/articles/s41586-024-07601-y)).       
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC ### `Protein` Language Model (`pLM`): Classify `Membrane Transport` Proteins
# MAGIC
# MAGIC For the purpose of our [AI/BI Dashboard](https://docs.databricks.com/en/dashboards/index.html) with [enabled Genie Space](https://docs.databricks.com/en/dashboards/index.html#enable-a-genie-space-from-your-dashboard) example, we will perform classification of `membrane transport proteins` on our [DLT](https://e2-demo-west.cloud.databricks.com/pipelines/c6a3e57b-1c44-476f-9e56-c8e05d9975f5/updates/2adf25f7-deb6-40fa-8f87-68efd0ca05a0?o=2556758628403379%3Fparent%3Dfolders%2F1625373258638091) [extracted protein sequences](/Workspace/dbdemos/HLS-ai-drug-discovery/01_DLT_ProteinData_processing). 
# MAGIC
# MAGIC We will leverage `protein language model` based on the NLP technique Bidirectional Encoder Representations from Transformers ([BERT](https://github.com/google-research/bert?tab=readme-ov-file#introduction)) adapted for proteins to learn contextual embeddings of individual amino acids within a protein sequence, i.e. transformer based `protein Language Models` pretrained on protein sequences ([`ProtTrans`](https://ieeexplore.ieee.org/document/9477085): [`ProtBERT-BFD`](https://huggingface.co/Rostlab/prot_bert_bfd)) 
# MAGIC
# MAGIC Specifically, we will use -- [**`Rostlab/prot_bert_bfd_membrane`**](https://huggingface.co/Rostlab/prot_bert_bfd_membrane/tree/main) -- which is `ProtBERT-BFD` fine-tuned on the membrane proteins dataset. 
# MAGIC
# MAGIC In a [recent research benchmarking a hybrid model approach](https://www.degruyter.com/document/doi/10.1515/jib-2022-0055/html) of combining [`ProtTrans`](https://ieeexplore.ieee.org/document/9477085) (e.g. [`Rostlab/prot_bert_bfd`](https://huggingface.co/Rostlab/prot_bert_bfd)) with convolutional neural networks (`CNN`), the **`prot_bert_bfd_membrane` + `CNN`** combination was shown to classify membrane transport protein using protein sequence with cross-validated sensitivity and specificity scores of ~95+% and independent test set sensitivity and specificity scores of ~91%.
# MAGIC
# MAGIC <!-- bfd: big fantastic dataset https://bfd.mmseqs.com/ -->

# COMMAND ----------

# MAGIC %md
# MAGIC ---    

# COMMAND ----------

# MAGIC %md
# MAGIC ## Protein Classification Setup & Batch Inference

# COMMAND ----------

# MAGIC %md
# MAGIC ####[1] Install and import relevant libraries .e.g `torch` and `transformers` 

# COMMAND ----------

# DBTITLE 1,Install torch & transformers
# MAGIC %pip install -q torch==2.3.1 transformers==4.41.2
# MAGIC
# MAGIC dbutils.library.restartPython()
# MAGIC

# COMMAND ----------

# DBTITLE 1,run nb utils
# MAGIC %run ./utils

# COMMAND ----------

# DBTITLE 1,Import torch / tokenizer / transformers pipelines
import torch 
from transformers import AutoTokenizer, AutoModelForSequenceClassification, TextClassificationPipeline
import re

# COMMAND ----------

# MAGIC %md
# MAGIC ####[2] Initialize the `prot_bert_bfd_membrane` transformer model hosted on Hugging Face
# MAGIC
# MAGIC - https://huggingface.co/Rostlab/prot_bert_bfd_membrane

# COMMAND ----------

# DBTITLE 1,Define TextClassificationPipeline
pipeline = TextClassificationPipeline(
    model=AutoModelForSequenceClassification.from_pretrained("Rostlab/prot_bert_bfd_membrane"),
    tokenizer=AutoTokenizer.from_pretrained("Rostlab/prot_bert_bfd_membrane"),
    device='cuda' if torch.cuda.is_available() else 'cpu' ## it is more efficient to run this with Classic GPU compute. 
)

# COMMAND ----------

# MAGIC %md
# MAGIC ####[3] Initial model inferencing before bulk inferencing
# MAGIC
# MAGIC - The model requires that there is a space between each capitalized amino acid
# MAGIC

# COMMAND ----------

# DBTITLE 1,input the spaced_sequence into the pipeline()
pipeline("M E S N L S G L V P A A G L V P A L P P A V T L G L T A A Y T T L Y A L L F F S V Y A Q L W L V L L Y G H K R L S Y Q T V F L A L C L L W A A L R T T L F S F Y F R D T P R A N R L G P L P F W L L Y C C P V C L Q F F T L T L M N L Y F A Q V V F K A K V K R R P E M S R G L L A V R G A F V G A S L L F L L V N V L C A V L S H R R R A Q P W A L L L V R V L V S D S L F V I C A L S L A A C L C L V A R R A P S T S I Y L E A K G T S V C Q A A A M G G A M V L L Y A S R A C Y N L T A L A L A P Q S R L D T F D Y D W Y N V S D Q A D L V N D L G N K G Y L V F G L I L F V W E L L P T T L L V G F F R V H R P P Q D L S T S H I L N G Q V F A S R S Y F F D R A G H C E D E G C S W E H S R G E S T R C Q D Q A A T T T V S T P P H R R D P P P S P T E Y P G P S P P H P R P L C Q V C L P L L A Q D P G G R G Y P L L W P A P C C S C H S E L V P S P")

# COMMAND ----------

# MAGIC %md
# MAGIC ####[4] Downsample for batch inferencing 
# MAGIC
# MAGIC For this demo, we work with a "_tiny_" sample of 500 proteins to do inference instead of the full 500,000 proteins. 
# MAGIC We can scale our solution for a production workload.

# COMMAND ----------

# DBTITLE 1,PLEASE UPATE
remove_widgets() 
uc_config = setup_uc_paths(spark=None, use_widgets=False); ## if you update the values in widgets -- it will automatically trigger an update of the UC paths

# Extract catalog, schema, volume names
catalog_name = uc_config["catalog_name"]
schema_name = uc_config["schema_name"]
volume_name = uc_config["volume_name"]

# COMMAND ----------

# DBTITLE 1,downsample full dataset
# Set the schema
spark.sql(f"USE {catalog_name}.{schema_name}")

# Create a sample dataset from large data for demo purposes ~ 500 records -- increase if desired
spark.sql(f"""
    CREATE OR REPLACE TABLE {catalog_name}.{schema_name}.tiny_sample_data AS
    SELECT *
    FROM {catalog_name}.{schema_name}.enriched_protein
    TABLESAMPLE (0.1 PERCENT)
""")

# COMMAND ----------

# MAGIC %md
# MAGIC ####[5] Vectorize Protein Classification 
# MAGIC We will demonstrate how we vectorize the protein classification by defining a [Pandas UDF](https://docs.databricks.com/en/udf/pandas.html) and applying it to our protein sequences as a batch process.

# COMMAND ----------

# DBTITLE 1,PandasUDF for Protein Classification
# import pyspark.sql class functions and types
from pyspark.sql import functions as F, types as T
import pandas as pd

schema = T.StructType([
    T.StructField("label", T.StringType(), True),
    T.StructField("score", T.FloatType(), True)
])

@F.pandas_udf(schema)
def classify_protein(sequences: pd.Series) -> pd.DataFrame:
    results = [pipeline(sequence) for sequence in sequences]
    labels = [result[0]['label'] for result in results]
    scores = [result[0]['score'] for result in results]
    return pd.DataFrame({
        "label": labels,
        "score": scores
    })

# COMMAND ----------

# DBTITLE 1,Apply UDF to Protein sparkDF
df = spark.read.table(f"{catalog_name}.{schema_name}.tiny_sample_data")

## The model requires a space between each amino acid -- we can add a space to the input sequence using spark functions `expr`        
df = df.withColumn('spaced_sequence', F.expr("concat_ws(' ', split(sequence, ''))"))
df = df.withColumn('spaced_sequence', F.expr("trim(spaced_sequence)"))


# Assuming df is your DataFrame with the 'spaced_sequence' column
# Apply the UDF to the 'spaced_sequence' column and create new columns for label and score
df = df.withColumn("classification", classify_protein("spaced_sequence"))

df = df.select("*", "classification.*") # extract the json struct into columns

display(df.limit(100))

# COMMAND ----------

# MAGIC %md
# MAGIC ####[6] Write out classified proteins sparkDF to Unity Catalog 
# MAGIC We will use the classified proteins in the downstream [AI/BI Dashboard](https://docs.databricks.com/en/dashboards/index.html) with [enabled Genie Space](https://docs.databricks.com/en/dashboards/index.html#enable-a-genie-space-from-your-dashboard) example implementation.

# COMMAND ----------

# DBTITLE 1,Write out the output to UC catalog
df.write.mode("overwrite").option("mergeSchema", "true").saveAsTable(f"{catalog_name}.{schema_name}.ProteinClassification_tiny")
