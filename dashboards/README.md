# AI-Driven Drug Discovery AIBI Genie 

### Dashboard Setup

The dashboard template JSON file:    
`./dashboards/AI-DrivenDrugDiscovery_[Overview&Dashboard<catalog_name>.ai_driven_drug_discovery].lvdash.json` includes the corresponding queries and code needed to generate a similar dashboard.

- Please update all Unity Catalog **`{catalog}.{schema}.*table`** paths to match the ones being specified and used in the notebooks, e.g. users may have specified their **`<catalog_name>.<schema_name>`** etc. Default `Table` names would be used per notebook/dashboard template. 

- Reivew and `run` the notebook [**`code2embed_img_in_dashboard.ipynb`**](code2embed_img_in_dashboard.ipynb) to copy over the relevant architecture flow image for embedding in dashboard, remembering to update the corresponding `markdown-url` in the dashboard template `JSON` as shown in the notebook.

- For "ask Genie" within dashboard, please **Enable Genie** via the **`AIBI settings`** (see screenshot below; RHS panel). Alternatively you can create a separate Genie Space and link it to the dashboard.      

- Additional [example questions](dashboards/GenieSpace_AdditionalExampleQs_GeneratedQueries.ipynb) are provided as guidance. 
  <br>

  ![AIBI settings](../assets/imgs/Enable_GenieSpace_via_AIBI_settings.png)

<br>

---    
<br>

### What the 2-page Dashboard looks like:

#### [Overview Tab] -- `print-extracted from editing view` 
![overview](../assets/imgs/AI-Driven%20Drug%20Discovery%20%5BOverview%5D.jpg)    
<!-- usingg github filename instead  -->
   
#### [Dashboard Tab] -- `downloaded as PDF from published view`
![dashboard](../assets/imgs/AI-Driven%20Drug%20Discovery%20%5BDashboard%5D.jpg)

