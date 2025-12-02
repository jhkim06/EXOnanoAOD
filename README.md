# EXOnanoAOD
Development of custom cms EXO PAG nanoAOD format


## Setup in CMSSW
```
cmsrel CMSSW_16_0_0_pre1
cd CMSSW_16_0_0_pre1/src
cmsenv
git cms-init
git cms-addpkg PhysicsTools/NanoAOD
cd PhysicsTools
git clone git@github.com:kerstinlovisa/EXOnanoAOD.git
cd ..
scram b -j
```

# Setup
Please include your customizations as python scripts under `python`, and if you need a custom producer add it under `plugins`.

We advice to use the templates in the branch `EXOnanoAOD_template` as base and create a new branch from there. 
```
git checkout EXOnanoAOD_template
git checkout -b your_branch_name
```
Rename the plugin in `plugins/EXOnanoAODProducerTemplate.cc` and python script in `python/custom_exonanoaod_template_cff.py`, and fill them in for your customizations. Edit the test script in `test/Run3_2023_PAT_EXONANO_template.py` to include your customizations under `# EXOnanoAOD customisation`. 

# Test run AOD -> EXOnanoAOD for Run 3
In `test` there is a test config file for Run 3 2023 which runs over one MC file (from TTto2L2Nu dataset). It will run AOD to NanoAOD format directly. 

After adding your customization functions at the end of the script under `# EXOnanoAOD customisation` you should be able to run your customizations with:
```
cmsRun Run3_2023_PAT_EXONANO_template.py
```

# Event size of customized NanoAOD
Check the event size of your custom EXOnanoAOD implementations in your nanoAOD root file (replace `nanoAOD.root`) by running:
```
git cms-addpkg PhysicsTools/NanoAOD
$CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py nanoAOD.root --size size.html
```
To check multiple datasets, edit use the script: 
```
sh test/run.sh
```
where output json files can be analyzed by `test/analyze.py`
