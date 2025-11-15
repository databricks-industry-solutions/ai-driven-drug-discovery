# AI-Driven Drug Discovery
How GenAI can help identify promising, under-researched proteins and uncover their properties for drug discovery

[![Databricks](https://img.shields.io/badge/Databricks-Solution_Accelerator-FF3621?style=for-the-badge&logo=databricks)](https://databricks.com)
[![Unity Catalog](https://img.shields.io/badge/Unity_Catalog-Enabled-00A1C9?style=for-the-badge)](https://docs.databricks.com/en/data-governance/unity-catalog/index.html)
[![Serverless](https://img.shields.io/badge/Serverless-Compute-00C851?style=for-the-badge)](https://docs.databricks.com/en/compute/serverless.html)

## Overview  

---   
### [todo update links]

The associated [set of notebooks](path?) corresponding to the [**_AI-Driven Drug Discovery_**](path?) mainly involves data preparation and  foundational model endpoint setup and usage. 

The notebooks cover: [01 protein data pre-processing](path?), [02 protein classification](path?), and will also leverage GenAI to [03 help democratize scientific knowledge](path?) as well as with [04 downstream discovery efforts](path?). 

The following diagram illustrates the flow of processes involved and the [Databricks Intelligent Platform](https://www.databricks.com/resources/demos/tours/horizontal/introducing-databricks-intelligence-platform) features utilized (in red).     

![architecture sketch](./assets/imgs/AI-Drug-Discovery-Page2.png) 

The [processed data written to Unity Catalog](https://docs.databricks.com/aws/en/data-governance/unity-catalog) and [served external foundational model endpoint](https://docs.databricks.com/en/generative-ai/external-models.html) will be used to create the associated [AI/BI Dashboard](https://docs.databricks.com/aws/en/dashboards/); e.g. [saved as pdf for sharing](./assets/imgs/`AI-Driven Drug Discovery [Dashboard].jpg`). 

<!-- > _**`to provide an example of how GenAI can be leveraged to discover more about protein properties and whether there are candidate proteins that are under-researched in areas that show promise for potential drug discovery.`**_ -->

**NB:** we provide an exported dashboard template in the `/dashboards` folder

When a [AIBI dashboard is Genie Space enabled](https://docs.databricks.com/aws/en/dashboards/genie-spaces) during publication, users can _"Ask Genie"_ about the datasets and [additional example questions](./dashboards/GenieSpace_AdditionalExampleQs_GeneratedQueries) are provided as guidance. 

<!-- Additional information on [resources used](https://e2-demo-west.cloud.databricks.com/editor/files/1625373258642799?o=2556758628403379) are included for reference. -->
   
<br> 
     
---    

Our example walkthrough highlights the following:    

#### Databricks Features

| Area | Feature | Role in this solution | Why it matters |
|------|---------|-----------------------|----------------|
| **Data & Governance** | **[Lakehouse / Delta](https://docs.databricks.com/en/delta/index.html)** | Stores FASTA raw files and derived protein tables (`enriched_protein`, semantic and research-enriched tables). | Central, scalable store for all protein data and features. |
| **Data & Governance** | **[Unity Catalog](https://docs.databricks.com/en/data-governance/unity-catalog/index.html)** | Governs access to tables and SQL functions (`scientific2simple`, `get_protein_research_info`). | Secure, governed foundation shared by DLT, AI/BI dashboards, and Genie. |
| **Data Engineering** | **[Delta Live Tables (DLT)](https://docs.databricks.com/en/delta-live-tables/index.html)** | Pipeline that ingests and preprocesses FASTA into structured protein tables with quality checks. | Reliable, declarative ETL so downstream AI/BI always sees clean, current data. |
| **Computation** | **[pandas UDFs](https://docs.databricks.com/en/sql/user-defined-functions/python.html)** | Applies vectorized protein language model inference during the Protein Classification step. | Scales Python/ML inference across large protein datasets. |
| **SQL & AIBI** | **[Databricks SQL](https://docs.databricks.com/en/sql/index.html)** | Runs aggregations and joins for distributions, counts, and high-confidence protein subsets. | Single query layer for both classic analytics and AI-enriched columns. |
| **SQL & AIBI** | **[AI/BI Dashboards](https://docs.databricks.com/en/ai-bi/index.html)** | Presents tiles for summary stats, distributions, organism comparisons, and detailed protein tables. | Visual exploration of protein data for non-coders. |
| **SQL & AIBI** | **Dashboard parameters & visuals** | Filters (Organism_SimpleTerm, score thresholds, ProteinType) and charts/tables bound to SQL. | Interactive “what‑if” exploration within the dashboard. |
| **AI Integration** | **[`ai_query()`](https://docs.databricks.com/aws/en/sql/language-manual/functions/ai_query) / [AI functions](https://docs.databricks.com/aws/en/large-language-models/ai-functions)** | SQL entry point to LLMs, used inside `scientific2simple` and `get_protein_research_info`. | Directly augments tables with AI-generated fields from within SQL. |
| **AI Integration** | **[Databricks hosted LLMs as FMAPI (e.g., Llama 3)](https://docs.databricks.com/aws/en/machine-learning/model-serving/score-foundation-models)** | Power `scientific2simple()` to convert scientific organism/protein text into layman terms. | Creates the semantic layer (`Organism_SimpleTerm`, `SimpleTerms`). |
| **AI Integration** | **[External models (e.g., GPT‑4o)](https://docs.databricks.com/en/generative-ai/external-models.html)** | Power `get_protein_research_info()` to attach recent/under‑researched context. | Turns rows into actionable research leads. |
| **AI Integration** | **[UC SQL Functions](https://docs.databricks.com/en/sql/language-manual/sql-ref-syntax-ddl-create-function.html)** | Registers `scientific2simple()` and `get_protein_research_info()` as reusable SQL functions. | Standardizes and reuses AI logic across queries, dashboards, and Genie. |
| **Exploration UX** | **EDA & ad‑hoc queries** | Parameterized tiles and queries for distributions, diversity, and high‑confidence subsets. | Guided and custom exploration on the same governed data. |
| **Exploration UX** | **Precomputed enrichment** | Stores AI-enriched research info in UC tables consumed by the dashboard. | Fast interaction without repeated external model calls. |
| **Conversational Analytics** | **[Genie](https://docs.databricks.com/en/genie/index.html)** | Natural-language interface over the same UC tables and functions. | Lets users explore proteins and insights via chat instead of SQL. |
| **Operationalization** | **AI semantic layer** | `Organism_SimpleTerm` and `SimpleTerms` reused by dashboards and Genie. | Bridges scientific jargon and business-friendly search/filter terms. |



<!-- ### Databricks Features

| Area | Feature | How it’s used in this Solution Accelerator | Why it matters / Highlight |
|------|---------|-----------------------------------------------|----------------------------|
| **Data & Governance** | **[Delta Lake (Delta tables)](https://docs.databricks.com/en/delta/index.html)** | FASTA protein information, sequences, classifications, and enriched outputs (e.g., `enriched_protein`, `organism_info_scientificnsimple`, `organism_protein_research_info`) are stored as Delta tables in the Lakehouse. | Provides a single source of truth for all protein data and derived features, enabling consistent BI + AI workloads on governed, scalable storage. |
| **Data & Governance** | **[Unity Catalog (UC)](https://docs.databricks.com/en/data-governance/unity-catalog/index.html)** | Governs access to raw FASTA files and all derived tables, plus AI-backed SQL functions like `scientific2simple()` and `get_protein_research_info()`. | Ensures secure, governed access to both data and AI functions so the same assets safely power AI/BI Dashboards and Genie. |
| **Data Engineering** | **[Delta Live Tables (DLT)](https://docs.databricks.com/en/delta-live-tables/index.html)** | A DLT pipeline ingests and preprocesses FASTA raw files into structured protein tables (e.g., `enriched_protein`), applying quality checks and transformations. | Provides a declarative, reliable pipeline for building and maintaining up-to-date protein datasets that the downstream AI and dashboards can trust. |
| **Computation / Transformation** | **[pandas UDFs](https://docs.databricks.com/en/sql/user-defined-functions/python.html)** | A pandas UDF is used in the Protein Classification step to apply vectorized protein language model inference over batches of sequences. | Enables scalable, high-throughput application of ML inference or complex Python logic directly over large data sets in the Lakehouse. |
| **SQL & AIBI** | **[Databricks SQL](https://docs.databricks.com/en/sql/index.html)** | Primary query engine for aggregations and EDA: computing distributions of molecular weight, counts by ProteinType and Organism, filtering by `ProteinClassificationScore`, and joining to AI‑generated columns. | Unifies traditional analytics and AI-enriched columns in a single SQL layer, making advanced analysis accessible to SQL users. |
| **SQL & AIBI** | **[AI/BI Dashboards](https://docs.databricks.com/en/ai-bi/index.html)** | The “HLS - AI Driven Drug Discovery” UI is built with AI/BI Dashboards, with tiles for summary statistics, EDA, comparative analyses, and detailed tables. | Delivers a modern BI experience tightly integrated with Databricks AI capabilities, so domain experts can visually explore insights without writing code. |
| **SQL & AIBI** | **[Dashboard parameters / controls](https://docs.databricks.com/en/ai-bi/index.html)** | Interactive controls allow selection of `Organism_SimpleTerm`, `ProteinClassificationScore` thresholds, and `ProteinType` (Soluble / Membrane). These are bound to SQL queries powering each tile. | Lets users ask “what if?” questions (e.g., raising confidence thresholds, changing organisms) in real time through simple UI controls in AI/BI Dashboards. |
| **SQL & AIBI** | **[Visualizations](https://docs.databricks.com/en/ai-bi/index.html)** | Uses histograms/box plots for molecular weight distributions, bar charts for protein-type percentages and organism frequencies, and tables for detailed protein attributes. | Makes complex protein data interpretable at a glance in the AI/BI Dashboard UX, surfacing patterns like type balance and weight distributions across organisms. |
| **AI & ML Integration** | **[`ai_query()`](https://docs.databricks.com/aws/en/sql/language-manual/functions/ai_query) / [AI functions](https://docs.databricks.com/aws/en/large-language-models/ai-functions)** | Called from SQL in functions like `scientific2simple()` and `get_protein_research_info()` to convert scientific descriptors to simple terms and to fetch research-related information. | Bridges Lakehouse data and LLMs directly from SQL, enabling AI-powered enrichment and transformations inside normal queries. |
| **AI & ML Integration** | **[Databricks-hosted LLMs as FMAPI (e.g., Llama 3)](https://docs.databricks.com/aws/en/machine-learning/model-serving/score-foundation-models)** | `scientific2simple()` uses `ai_query()` with a Databricks-hosted Llama 3 model to generate `Organism_SimpleTerm` and `SimpleTerms` arrays. | Creates an AI-driven semantic layer that turns jargon-heavy organism/protein names into terms non-experts can search and filter by. |
| **AI & ML Integration** | **[External models (e.g., Azure OpenAI GPT‑4o)](https://docs.databricks.com/en/generative-ai/external-models.html)** | `get_protein_research_info()` calls `ai_query()` against an external GPT‑4o endpoint to enrich proteins with recent_research and under_researched-area information. | Transforms the dashboard from a static catalog into a research tool that highlights where scientific opportunities may exist. |
| **AI & ML Integration** | **[UC-registered SQL functions](https://docs.databricks.com/en/sql/language-manual/sql-ref-syntax-ddl-create-function.html)** | Functions like `scientific2simple()` and `get_protein_research_info()` are registered in UC and used directly in SQL queries and dashboards. | Operationalizes LLM usage as standard SQL, making AI enrichment consistent, auditable, and easy to reuse across teams and workloads. |
| **Exploratory & Analytical UX** | **Interactive EDA in SQL / AI/BI dashboards** | Sections “Summary Statistics”, “EDA”, and “Protein Insights” are implemented as parameterized tiles showing distributions, high-confidence subsets, organism diversity, and shared proteins. | Provides a structured exploration path: from “what do we have?” to “where are high-confidence, potentially interesting protein targets?”. |
| **Exploratory & Analytical UX** | **Ad‑hoc parameterized queries** | Section 3.1.1 supports ad-hoc, parameterized queries using `Organism_SimpleTerm` and `ProteinClassificationScore`, including more expensive AI-enriched queries. | Gives power users a way to go deeper than canned tiles, tailoring high-cost AI analyses to specific organisms and thresholds. |
| **Exploratory & Analytical UX** | **Precomputation for performance** | Section 3.1.2 uses precomputed research info (e.g., `organism_protein_research_info`) surfaced via `get_protein_research_info()` for fast exploration. | Balances insight and performance: users get AI-enriched context quickly without repeated live calls to external models. |
| **Conversational Analytics** | **[Genie (AI/BI Genie)](https://docs.databricks.com/en/genie/index.html)** | A Genie space is associated with the same UC data and AI functions; users ask natural-language questions and Genie generates SQL, tables, or charts. | Enables conversational analytics on top of the same Lakehouse — non-SQL users can still explore protein data and insights using natural language. |
| **Operationalization & Reuse** | **AI-driven semantic layer (Organism_SimpleTerm & SimpleTerms)** | AI-generated simple terms (via `scientific2simple()`) represent organisms and concepts in layman language and are reused by both dashboards and Genie. | Creates a reusable semantic layer that connects scientific data with business-friendly language, driving both UX and discoverability. | -->

---   

## Getting Started

Clone this repository to your Databricks Workspace.  
You will find the set of notebooks referenced above within the `/notebooks` folder, the dashboards json files within `/dashboards` -- UC-tables will need to be re-referenced after they are being generated via the notebooks.  

---     

## Contributing

**We welcome contributions!**  

1. **git clone** this project locally
2. Utilize the Databricks CLI to test your changes against a Databricks workspace of your choice
3. Contribute to repositories with pull requests (PRs), ensuring that you always have a second-party review from a capable teammate

<!-- Please refer to [REPO_structure.md](REPO_structure.md) and [CONTRIBUTING.md](CONTRIBUTING.md) for more details and guidance.     -->

---   

## How to get help
Databricks support doesn't cover this content. For questions or bugs, please open a GitHub issue and the team will help on a best effort basis.   

---   

## Licenses

&copy; 2025 Databricks, Inc. All rights reserved. The source in this project is provided subject to the Databricks License [https://databricks.com/db-license-source]. All included or referenced third party libraries are subject to the licenses set forth below.

| Package | License | Copyright |
|---------|---------|-----------|
| [biopython](https://github.com/biopython/biopython) | [Biopython License Agreement (dual licensed with BSD 3-Clause)](https://github.com/biopython/biopython/blob/master/LICENSE.rst) | Copyright (c) 1999-2024, The Biopython Contributors |
| [databricks-sdk](https://github.com/databricks/databricks-sdk-py) | [Apache License 2.0](https://github.com/databricks/databricks-sdk-py/blob/main/LICENSE) | Copyright 2023-present, Databricks, Inc. |
| [pandas](https://github.com/pandas-dev/pandas) | [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/main/LICENSE) | Copyright (c) 2008-2011, AQR Capital Management, LLC, Lambda Foundry, Inc. and PyData Development Team; Copyright (c) 2011-present, Open source contributors |
| [pyspark](https://github.com/apache/spark) | [Apache License 2.0](https://github.com/apache/spark/blob/master/LICENSE) | Copyright 2014-present, The Apache Software Foundation |
| [Rostlab/prot_bert_bfd_membrane](https://huggingface.co/Rostlab/prot_bert_bfd_membrane) | [Academic Free License v3.0 (Model Weights)](https://opensource.org/license/afl-3-0-php) | Copyright 2020, Elnaggar et al. (Rostlab, TUM) |
| [torch](https://github.com/pytorch/pytorch) | [BSD 3-Clause License (Modified BSD)](https://github.com/pytorch/pytorch/blob/main/LICENSE) | Copyright (c) 2016-present, Facebook, Inc. and contributors |
| [transformers](https://github.com/huggingface/transformers) | [Apache License 2.0](https://github.com/huggingface/transformers/blob/main/LICENSE) | Copyright 2018-present, The Hugging Face team |
| [UniProt (uniprot_sprot.fasta.gz)](https://www.uniprot.org/help/downloads) | [Creative Commons Attribution 4.0 (CC BY 4.0)](https://www.uniprot.org/help/license) | Copyright 2002-present, The UniProt Consortium |


<!-- | Package | License | Copyright |
|---------|---------|-----------|
| | | | -->

---   
