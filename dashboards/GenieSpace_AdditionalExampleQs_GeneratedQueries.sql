-- Databricks notebook source
-- MAGIC %md
-- MAGIC ### Here are some additional `Example Queries` you can try asking Genie... 

-- COMMAND ----------

-- MAGIC %md
-- MAGIC Chat Examples 
-- MAGIC - Research potential of protein types 
-- MAGIC   - [For each of the top 10 membrane transporter and soluble proteins look up their research potential and output with simple terms](https://e2-demo-west.cloud.databricks.com/genie/rooms/01efc9fc59e2146eab09fc4e6a4b9a5a/chats/01efd78c0b5217c489b0836cc9ca2b33?m=01efd78c8aa0119c874f80ff0af2bbb9) 
-- MAGIC
-- MAGIC - Is it possible to search for common proteins with similar under researched potential?
-- MAGIC   - [Search for common proteins with similar under researched potential for cancer drugs and associated organism simple terms; remove duplicates](https://e2-demo-west.cloud.databricks.com/genie/rooms/01efc9fc59e2146eab09fc4e6a4b9a5a/chats/01efd864563913178d5a2a7ff75d83e6?m=01efd86475c41537b3a7204a3b89c3de)

-- COMMAND ----------

-- DBTITLE 1,common proteins shared by multiple organism
-- # - Are there common proteins shared by multiple organisms? Can you summarize with a count of a number of organisms sharing common proteins including their simple names and a corresponding list of Organism names and the common Proteins shared along with the protein information and its classification score?

-- Generated Query
WITH protein_organism AS (
  SELECT
    `ProteinName`,
    `OrganismName`,
    `score`
  FROM
    `demos_genie`.`hls_ai_drug_discovery`.`proteinclassification_tiny`
  WHERE
    `ProteinName` IS NOT NULL
    AND `OrganismName` IS NOT NULL
    AND `score` IS NOT NULL
),
protein_count AS (
  SELECT
    `ProteinName`,
    COUNT(DISTINCT `OrganismName`) AS `Organism_Count`
  FROM
    protein_organism
  GROUP BY
    `ProteinName`
  HAVING
    COUNT(DISTINCT `OrganismName`) > 1
)
SELECT
  p.`ProteinName`,
  p.`OrganismName`,
  p.`score`,
  o.`Organism_SimpleTerm`,
  o.`Organism_Definition`
FROM
  protein_organism p
  JOIN protein_count pc ON p.`ProteinName` = pc.`ProteinName`
  JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_info_scientificnsimple` o ON p.`OrganismName` = o.`OrganismName`
ORDER BY
  p.`ProteinName`,
  p.`OrganismName`

-- COMMAND ----------

-- DBTITLE 1,recent research and under-researched-areas from these common proteins
-- # -- Look up recent research and under-researched-areas from these common proteins?

-- Generated Query
WITH protein_organism AS (
  SELECT
    `ProteinName`,
    `OrganismName`,
    `score`
  FROM
    `demos_genie`.`hls_ai_drug_discovery`.`proteinclassification_tiny`
  WHERE
    `ProteinName` IS NOT NULL
    AND `OrganismName` IS NOT NULL
    AND `score` IS NOT NULL
),
protein_count AS (
  SELECT
    `ProteinName`,
    COUNT(DISTINCT `OrganismName`) AS `Organism_Count`
  FROM
    protein_organism
  GROUP BY
    `ProteinName`
  HAVING
    COUNT(DISTINCT `OrganismName`) > 1
)
SELECT
  p.`ProteinName`,
  p.`OrganismName`,
  p.`score`,
  r.`recent_research`,
  r.`under_researched_areas`
FROM
  protein_organism p
  JOIN protein_count pc ON p.`ProteinName` = pc.`ProteinName`
  JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_protein_research_info` r ON p.`OrganismName` = r.`OrganismName`
  AND p.`ProteinName` = r.`ProteinName`
ORDER BY
  p.`ProteinName`,
  p.`OrganismName`

-- COMMAND ----------

-- DBTITLE 1,common organisms for protein types
-- # -- For the types of proteins classified, provide the list of most common organism in their simple terms; provide info. on protein info., types and their scores for the organisms

-- Generated Query
WITH q AS (
  SELECT
    p.OrganismName,
    o.Organism_SimpleTerm,
    p.ProteinName,
    p.label AS ProteinType,
    p.score AS ProteinClassificationScore
  FROM
    `demos_genie`.`hls_ai_drug_discovery`.`proteinclassification_tiny` p
    JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_info_scientificnsimple` o ON p.OrganismName = o.OrganismName
  WHERE
    p.ProteinName IS NOT NULL
    AND p.label IS NOT NULL
    AND p.score IS NOT NULL
)
SELECT
  `Organism_SimpleTerm`,
  `ProteinName`,
  `ProteinType`,
  `ProteinClassificationScore`
FROM
  q
ORDER BY
  `Organism_SimpleTerm`,
  `ProteinType`,
  `ProteinClassificationScore` DESC

-- COMMAND ----------

-- DBTITLE 1,common organisms for top 5 protein types
-- # -- For the types of proteins classified, provide the list of most common organisms in their simple terms; provide info. on protein info., types, and their scores for the organisms. Then visualize the top 5 most frequently occurring organisms in Membrane and Soluble Protein Types separately, sorted by order of descending frequency — provide info on the type of protein — the output needs to be for both Membrane and Soluble Protein types of proteins; provide 2 visualizations 

--  It cannot generate 2 viz so we break it into 2 Qs


-- For the types of proteins classified, provide the list of most common organisms in their simple terms; provide info. on protein info. , types, and their scores for the organisms. Then visualize the top 5 most frequently occurring organisms in Membrane Protein Types, sorted by order of descending frequency — provide info on the type of protein 

-- Generated Query
WITH q AS (
  SELECT
    p.OrganismName,
    o.Organism_SimpleTerm,
    p.ProteinName,
    p.label AS ProteinType,
    p.score AS ProteinClassificationScore
  FROM
    `demos_genie`.`hls_ai_drug_discovery`.`proteinclassification_tiny` p
    JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_info_scientificnsimple` o ON p.OrganismName = o.OrganismName
  WHERE
    p.ProteinName IS NOT NULL
    AND p.label IS NOT NULL
    AND p.score IS NOT NULL
    AND p.label ILIKE '%membrane%'
)
SELECT
  `Organism_SimpleTerm`,
  `ProteinName`,
  `ProteinType`,
  `ProteinClassificationScore`,
  COUNT(*) AS Frequency
FROM
  q
GROUP BY
  `Organism_SimpleTerm`,
  `ProteinName`,
  `ProteinType`,
  `ProteinClassificationScore`
ORDER BY
  Frequency DESC
LIMIT
  5


-- do similar for soluble proteins type

-- Generated Query
WITH q AS (
  SELECT
    p.OrganismName,
    o.Organism_SimpleTerm,
    p.ProteinName,
    p.label AS ProteinType,
    p.score AS ProteinClassificationScore
  FROM
    `demos_genie`.`hls_ai_drug_discovery`.`proteinclassification_tiny` p
    JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_info_scientificnsimple` o ON p.OrganismName = o.OrganismName
  WHERE
    p.ProteinName IS NOT NULL
    AND p.label IS NOT NULL
    AND p.score IS NOT NULL
)
SELECT
  `Organism_SimpleTerm`,
  `ProteinName`,
  `ProteinType`,
  `ProteinClassificationScore`
FROM
  q
WHERE
  `ProteinType` ILIKE '%soluble%'
ORDER BY
  `Organism_SimpleTerm`,
  `ProteinType`,
  `ProteinClassificationScore` DESC
LIMIT
  5

-- COMMAND ----------

-- DBTITLE 1,proteins worthy of further research focus
-- # -- Based on proteins with high classification score, recent research, under-researched areas, and potential promise, which candidate proteins and corresponding organisms, in simple terms, might be worth exploring for further research and development? Sort by protein type, protein name 

-- Generated Query
WITH high_confidence_proteins AS (
  SELECT
    p.ProteinName,
    p.OrganismName,
    o.Organism_SimpleTerm,
    p.label,
    p.score,
    r.recent_research,
    r.under_researched_areas
  FROM
    `demos_genie`.`hls_ai_drug_discovery`.`proteinclassification_tiny` p
    JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_protein_research_info` r ON p.ProteinName = r.ProteinName
    AND p.OrganismName = r.OrganismName
    JOIN `demos_genie`.`hls_ai_drug_discovery`.`tinysample_organism_info_scientificnsimple` o ON p.OrganismName = o.OrganismName
  WHERE
    p.score > 0.85
)
SELECT
  ProteinName,
  Organism_SimpleTerm,
  label,
  recent_research,
  under_researched_areas
FROM
  high_confidence_proteins
ORDER BY
  label,
  ProteinName
